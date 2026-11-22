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
  # The basic A₀/A₁ Beyn method cannot handle this; it needs higher moments.
end
