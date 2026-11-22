# Beyn.jl

Beyn.jl solves nonlinear eigenvalue problems

```math
T(\lambda)\,x = 0, \qquad x \neq 0,
```

where ``T`` is an ``n \times n`` matrix-valued function that is analytic inside
an ellipse (by default the unit circle). It uses Beyn's contour-integral method
to find every eigenvalue ``\lambda`` inside the ellipse, together with an
eigenvector. A rectangle can be given instead; the solver then searches the
least-area ellipse around it.

## Quick start

```julia
using Beyn, LinearAlgebra, Random

Random.seed!(1)
Q = Matrix(qr(randn(ComplexF64, 5, 5)).Q)
D(z) = Diagonal([z - 0.3, z + 0.5im, exp(z) - 1.5, z - 2, z^2 - 4])
T(z) = Q * D(z) * Q'

sol = Beyn.solver(T, 5)
sol.eigenvalues      # ≈ 0.3, -0.5i and log(1.5)
sol.eigenvectors     # one column per eigenvalue
```

Only three of the five diagonal entries have a root inside the unit circle,
so the solver returns three eigenvalues. Because `D` is rotated by the unitary
matrix `Q`, the eigenvectors are the matching columns of `Q`, up to scaling.

## Checking a result

Compute the residual for each eigenpair:

```julia
for (λ, x) in zip(sol.eigenvalues, eachcol(sol.eigenvectors))
    @show norm(T(λ) * x) / norm(x)
end
```

The singular values in `sol.subproblemsingularvalues` show how the rank was
chosen. There should be a clear gap between the values that were kept (one
per eigenvalue) and the rest, which sit at round-off level.

## Next steps

- [Algorithm](@ref) explains what the solver computes.
- [API](@ref) lists the arguments, options and return values.
