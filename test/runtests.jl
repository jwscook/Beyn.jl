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

  @testset "argument principle count" begin
    @test Beyn.countzeros(z -> Diagonal([exp(z) - 1.5, z^3 - 0.001, z - 2])) == 4
    @test Beyn.countzeros(z -> fill(complex(z)^20 - 0.5^20, 1, 1)) == 20
    @test Beyn.countzeros(z -> Diagonal([z - 2, z + 3])) == 0
  end
end
