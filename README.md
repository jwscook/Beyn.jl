# Beyn.jl

Solve nonlinear eigenvalue problems with **Beyn's contour-integral method**.

Given a matrix-valued function `T(z)` that is analytic near the unit circle,
`Beyn.solver` finds every eigenvalue `λ` inside the unit disc, that is, every
value where `T(λ)` is singular, together with an eigenvector `x` satisfying
`T(λ) x = 0`.

## Installation

The package is not registered, so install it from its repository URL or a local path:

```julia
using Pkg
Pkg.add(url = "https://github.com/OWNER/Beyn.jl")   # or: Pkg.develop(path = "path/to/Beyn.jl")
```

Requires Julia 1.12 or newer.

## Quick start

```julia
using Beyn, LinearAlgebra, Random

Random.seed!(1)
Q = Matrix(qr(randn(ComplexF64, 5, 5)).Q)
D(z) = Diagonal([z - 0.3, z + 0.5im, exp(z) - 1.5, z - 2, z^2 - 4])
T(z) = Q * D(z) * Q'            # eigenvalues inside |z| < 1: 0.3, -0.5i, log(1.5)

sol = Beyn.solver(T, 5)          # use 5 random probe vectors

sol.eigenvalues                  # ≈ [0.3, -0.5im, log(1.5)] in some order
sol.eigenvectors                 # one column per eigenvalue
norm(T(sol.eigenvalues[1]) * sol.eigenvectors[:, 1])   # ≈ 0
```

## API

```julia
Beyn.solver(T, nprobes::Integer; kwargs...)
Beyn.solver(T, probes::AbstractMatrix; kwargs...)
```

- `T(z)` returns an `n × n` matrix for complex `z`. **It must return a new
  matrix on every call**, because the solver factorizes it in place.
- `nprobes` is the number of random, normalized probe vectors to use. Pass
  your own `n × ℓ` matrix `probes` instead to control them.

Keyword arguments:

| keyword    | default      | meaning                                                         |
|------------|--------------|-----------------------------------------------------------------|
| `Kmax`     | `4`          | max block rows of the Hankel matrices (stores `2Kmax` moments)  |
| `quadrtol` | `sqrt(eps())`| relative tolerance of the contour quadrature                    |
| `quadatol` | `0.0`        | absolute tolerance of the contour quadrature                    |
| `svdrtol`  | `sqrt(eps())`| singular values below `svdrtol * σ₁` are treated as zero        |
| `svdatol`  | `0.0`        | absolute cutoff for singular values                             |

The solver returns a named tuple:

| field                       | contents                                                  |
|-----------------------------|-----------------------------------------------------------|
| `eigenvalues`               | eigenvalues found inside the unit circle                  |
| `eigenvectors`              | matrix whose columns are the matching eigenvectors        |
| `subproblemsingularvalues`  | singular values of `H₀`, useful for checking the rank cut |
| `count`                     | eigenvalues inside (with multiplicity), by argument principle |
| `Keig`                      | block rows of the Hankel matrices used, `≤ Kmax`          |

## Limitations

- **Contour:** only the unit circle is supported. To search another disc
  `|z - c| < ρ`, pass `z -> T(c + ρ * z)` and map the results back with `c + ρ * λ`.
- **Entire `T`:** the eigenvalue count `m` is the winding number of `det T`
  around the contour, so `T` must have no poles inside it.
- **Number of eigenvalues:** with `ℓ` probes, the solver needs
  `Keig = cld(m, min(n, ℓ)) ≤ Kmax` block rows. It errors otherwise; raise
  `Kmax` or the number of probes. Block Hankel matrices also handle Jordan
  blocks, but `ℓ` must be at least the largest geometric multiplicity.
- **Eigenvalues near the contour:** accuracy drops as an eigenvalue gets close
  to `|z| = 1`, because the trapezoid rule then converges slowly.

## Documentation

Build the documentation locally with:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

## Reference

W.-J. Beyn, *An integral method for solving nonlinear eigenvalue problems*,
Linear Algebra and its Applications 436 (2012) 3839–3863.
