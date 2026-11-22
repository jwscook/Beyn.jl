# API

## `Beyn.solver`

```julia
Beyn.solver(T, nprobes::Integer; kwargs...)
Beyn.solver(T, probes::AbstractMatrix; kwargs...)
```

Find the eigenvalues of ``T(\lambda)x = 0`` inside the unit circle.

**Arguments**

- `T`: a function that takes a complex number `z` and returns an `n × n`
  matrix. It must return a new matrix on every call, because the solver
  overwrites it with its LU factorization.
- `nprobes`: the number of random probe vectors. They are drawn from a
  complex normal distribution and normalized. Use at least as many as the
  number of eigenvalues you expect inside the circle.
- `probes`: your own `n × ℓ` probe matrix, used instead of random probes.

**Keyword arguments**

| keyword    | default       | meaning                                                  |
|------------|---------------|----------------------------------------------------------|
| `quadrtol` | `sqrt(eps())` | relative tolerance of the contour quadrature             |
| `quadatol` | `0.0`         | absolute tolerance of the contour quadrature             |
| `svdrtol`  | `sqrt(eps())` | relative cutoff for the singular values of ``A_0``       |
| `svdatol`  | `0.0`         | absolute cutoff for the singular values of ``A_0``       |

**Returns** a named tuple with:

| field                      | contents                                                       |
|----------------------------|----------------------------------------------------------------|
| `eigenvalues`              | vector of eigenvalues inside the unit circle                   |
| `eigenvectors`             | matrix whose `i`th column is the eigenvector for `eigenvalues[i]` |
| `subproblemsingularvalues` | all singular values of ``A_0``, largest first                  |

The eigenvectors are not normalized.

**Errors**

- `"Maximum quadrature points reached"` if the quadrature has not converged
  after 65536 points. This usually means an eigenvalue lies on or very close
  to the unit circle.
- `SingularException` if `T(z)` is exactly singular at a quadrature point.

## Searching a different disc

Only the unit circle is built in. To search the disc ``|z - c| < \rho``,
rescale the problem and map the eigenvalues back:

```julia
c, ρ = 2.0, 0.5
sol = Beyn.solver(z -> T(c + ρ * z), 5)
λ = c .+ ρ .* sol.eigenvalues
```

The eigenvectors need no change.
