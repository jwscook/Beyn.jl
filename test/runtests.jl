using Beyn, Test, LinearAlgebra, Random

# f(λ)x ≈ 0, relative to ‖x‖
residual(f, λ, x) = norm(f(λ) * x) / norm(x)

@testset "Beyn" begin
  @testset "rotated diagonal nonlinear problem" begin
    Random.seed!(1)
    Q = Matrix(qr(randn(ComplexF64, 5, 5)).Q)
    D(z) = Diagonal([z - 0.3,          # 0.3       inside
                     z + 0.5im,        # -0.5i     inside
                     exp(z) - 1.5,     # log(1.5)  inside (log(1.5)+2πik, k≠0, outside)
                     z - 2,            # 2         outside
                     z^2 - 4])         # ±2        outside
    f(z) = Q * D(z) * Q'

    sol = Beyn.solver(f, 5)
    λ, X = sol.eigenvalues, sol.eigenvectors

    expected = [0.3, -0.5im, log(1.5)]
    @test length(λ) == 3
    @test all(any(isapprox(e, l; atol=1e-8) for l in λ) for e in expected)
    for i in eachindex(λ)
      @test residual(f, λ[i], X[:, i]) < 1e-8
    end
  end

  @testset "linear problem matches eigvals" begin
    Random.seed!(2)
    A = randn(ComplexF64, 6, 6) / 2
    f(z) = A - z * I
    sol = Beyn.solver(f, 6)
    inside = filter(z -> abs(z) < 1, eigvals(A))
    @test length(sol.eigenvalues) == length(inside)
    @test all(any(isapprox(e, l; atol=1e-8) for l in sol.eigenvalues) for e in inside)
  end

  @testset "no eigenvalues inside contour" begin
    f(z) = Diagonal([z - 2, z + 3])
    sol = Beyn.solver(f, 2)
    @test isempty(sol.eigenvalues)
  end

  # f(ω) = [ω²-1 2ω; -2ω ω²-1] has det = (ω²+1)²: each root ±i is a Jordan
  # block (algebraic multiplicity 2, geometric 1), total multiplicity 4 > n = 2.
  # The A₀/A₁ method cannot handle this; Keig = 2 block Hankel matrices can.
  # Rescaled by 2 so the roots ±i/2 are inside the unit circle.
  # Keig = 2 makes H₀ 4×4 with rank 4: saturated, so the solver grows to
  # Keig = 3 (6×6, rank 4), where a spare noise-level singular value confirms it.
  @testset "Jordan blocks, saturation grows Keig" begin
    Random.seed!(3)
    g(ω) = [ω^2-1 2ω; -2ω ω^2-1]
    f(z) = g(2z)
    sol = Beyn.solver(f, 2)
    @test sol.count == 4
    @test sol.Keig == 3
    @test sol.subproblemsingularvalues[5] < 1e-6 * sol.subproblemsingularvalues[1]
    @test length(sol.eigenvalues) == 4
    # a double root is only determined to ≈ √(perturbation)
    @test count(l -> isapprox(l, 0.5im; atol=1e-6), sol.eigenvalues) == 2
    @test count(l -> isapprox(l, -0.5im; atol=1e-6), sol.eigenvalues) == 2
    @test_throws ErrorException Beyn.solver(f, 2; Kmax=1)
    # Kmax = 2 still solves, but cannot confirm the count past saturation
    sol2 = @test_logs (:warn, r"saturated") Beyn.solver(f, 2; Kmax=2)
    @test sol2.Keig == 2
    @test count(l -> isapprox(l, 0.5im; atol=1e-6), sol2.eigenvalues) == 2
  end

  @testset "more eigenvalues than probes" begin
    Random.seed!(4)
    n = 8
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    vals = [0.1, -0.4, 0.5im, -0.3 - 0.3im, 0.6 + 0.2im, 1.5, -2, 3im]
    A = Q * Diagonal(vals) * Q'
    f(z) = A - z * I
    sol = Beyn.solver(f, 2)
    inside = filter(z -> abs(z) < 1, vals)
    @test sol.count == 5
    @test sol.Keig == 3
    @test length(sol.eigenvalues) == 5
    @test all(any(isapprox(e, l; atol=1e-8) for l in sol.eigenvalues) for e in inside)
    for i in eachindex(sol.eigenvalues)
      @test residual(f, sol.eigenvalues[i], sol.eigenvectors[:, i]) < 1e-8
    end
  end

  @testset "unit circle is the default ellipse" begin
    @test Beyn.Ellipse() == Beyn.Ellipse(0, 1, 1, 0)
    Random.seed!(5)
    A = randn(ComplexF64, 4, 4) / 2
    f(z) = A - z * I
    probes = randn(ComplexF64, 4, 4)
    s1 = Beyn.solver(f, probes)
    s2 = Beyn.solver(f, probes; contour=Beyn.Ellipse(0, 1))
    @test s1.eigenvalues ≈ s2.eigenvalues
    @test s1.contour == Beyn.Ellipse()
  end

  @testset "shifted, rotated, elongated ellipse" begin
    Random.seed!(6)
    n = 10
    e = Beyn.Ellipse(2 + 1im, 1.5, 0.3; angle=π / 6)
    # points at fractions of the way along each axis: inside if (x/a)²+(y/b)² < 1
    along(x, y) = Beyn.denormalize(e, complex(x * 1.5, y * 0.3) / 1.5)
    vals = [along(0, 0), along(0.8, 0), along(-0.5, 0.5), along(0.3, -0.7),   # inside
            along(0, 1.5), along(1.2, 0), along(-0.9, 0.6), 2 + 1im + 2im,    # outside
            along(0.75, 0.75), -1.0]                                          # outside
    @test count(in(e), vals) == 4
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    A = Q * Diagonal(vals) * Q'
    f(z) = A - z * I
    sol = Beyn.solver(f, 3; contour=e)
    inside = filter(in(e), vals)
    @test sol.count == 4
    @test length(sol.eigenvalues) == 4
    @test all(any(isapprox(v, l; atol=1e-8) for l in sol.eigenvalues) for v in inside)
    for i in eachindex(sol.eigenvalues)
      @test residual(f, sol.eigenvalues[i], sol.eigenvectors[:, i]) < 1e-8
    end
    @test Beyn.countzeros(f; contour=e) == 4
  end

  @testset "nonlinear problem in a shifted disc" begin
    Random.seed!(7)
    Q = Matrix(qr(randn(ComplexF64, 4, 4)).Q)
    # roots near 3: 3.1, 3 - 0.2i, log(20) ≈ 2.9957; 3.6 outside radius 0.5
    D(z) = Diagonal([z - 3.1, z - 3 + 0.2im, exp(z) - 20, z - 3.6])
    f(z) = Q * D(z) * Q'
    sol = Beyn.solver(f, 4; contour=Beyn.Ellipse(3, 0.5))
    expected = [3.1, 3 - 0.2im, log(20)]
    @test length(sol.eigenvalues) == 3
    @test all(any(isapprox(v, l; atol=1e-8) for l in sol.eigenvalues) for v in expected)
  end

  @testset "rectangle is enclosed by an ellipse" begin
    r = Beyn.Rectangle(1.2 + 0.2im, 0.5 - 0.1im)   # any two opposite corners
    @test r.lo == 0.5 - 0.1im && r.hi == 1.2 + 0.2im
    e = Beyn.Ellipse(r)
    @test e.centre ≈ 0.85 + 0.05im
    @test e.a ≈ √2 * 0.35 && e.b ≈ √2 * 0.15 && e.angle == 0
    # corners on the ellipse, edge midpoints inside
    w = cis(-e.angle) .* ([r.lo, r.hi, complex(real(r.lo), imag(r.hi))] .- e.centre)
    @test all(z -> (real(z) / e.a)^2 + (imag(z) / e.b)^2 ≈ 1, w)
    @test all(in(e), [0.5 + 0.05im, 1.2 + 0.05im, 0.85 - 0.1im, 0.85 + 0.2im])
    @test Beyn.Ellipse(r; margin=0.1).a ≈ 1.1 * e.a

    Random.seed!(8)
    n = 8
    vals = [0.6 + 0.0im, 1.0 + 0.1im, 0.9 - 0.05im,   # inside rectangle
            0.85 + 0.25im,                            # inside ellipse, outside rectangle
            0.0, 2.0, 0.85 + 0.5im, 1.0im]            # outside ellipse
    Q = Matrix(qr(randn(ComplexF64, n, n)).Q)
    A = Q * Diagonal(vals) * Q'
    f(z) = A - z * I
    sol = Beyn.solver(f, 4; contour=r)
    @test sol.contour == e
    @test sol.count == 4
    @test all(any(isapprox(v, l; atol=1e-8) for l in sol.eigenvalues) for v in vals[1:4])
    @test count(in(r), sol.eigenvalues) == 3
    @test Beyn.countzeros(f; contour=r) == 4
  end

  @testset "invalid contours" begin
    @test_throws ArgumentError Beyn.Ellipse(0, 0, 1)
    @test_throws ArgumentError Beyn.Ellipse(0, 1, -1)
    @test_throws ArgumentError Beyn.Ellipse(0, Inf)
    @test_throws ArgumentError Beyn.Rectangle(0, 1)        # zero height
    @test_throws ArgumentError Beyn.Rectangle(1im, 2im)    # zero width
    @test_throws ArgumentError Beyn.Ellipse(Beyn.Rectangle(0, 1 + 1im); margin=-0.1)
  end

  @testset "argument principle count" begin
    @test Beyn.countzeros(z -> Diagonal([exp(z) - 1.5, z^3 - 0.001, z - 2])) == 4
    @test Beyn.countzeros(z -> fill(complex(z)^20 - 0.5^20, 1, 1)) == 20
    @test Beyn.countzeros(z -> Diagonal([z - 2, z + 3])) == 0
  end
end
