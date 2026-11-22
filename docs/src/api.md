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
  complex normal distribution and normalized. Fewer probes than eigenvalues
  is fine as long as `Kmax` block rows make up the difference.
- `probes`: your own `n × ℓ` probe matrix, used instead of random probes.

**Keyword arguments**

| keyword    | default       | meaning                                                  |
|------------|---------------|----------------------------------------------------------|
| `Kmax`     | `4`           | maximum block rows of the Hankel matrices; moments ``A_0 \dots A_{2K_{\max}-1}`` are stored |
| `quadrtol` | `sqrt(eps())` | relative tolerance of the contour quadrature             |
| `quadatol` | `0.0`         | absolute tolerance of the contour quadrature             |
| `svdrtol`  | `sqrt(eps())` | relative cutoff for the singular values of ``H_0``       |
| `svdatol`  | `0.0`         | absolute cutoff for the singular values of ``H_0``       |

**Returns** a named tuple with:

| field                      | contents                                                       |
|----------------------------|----------------------------------------------------------------|
| `eigenvalues`              | vector of eigenvalues inside the unit circle                   |
| `eigenvectors`             | matrix whose `i`th column is the eigenvector for `eigenvalues[i]` |
| `subproblemsingularvalues` | all singular values of ``H_0``, largest first                  |
| `count`                    | eigenvalues inside, with multiplicity, from the argument principle |
| `Keig`                     | block rows used, ``\le K_{\max}`` (0 if `count == 0`)           |

The eigenvectors are not normalized. For an eigenvalue with a Jordan block,
the repeated copies of ``\lambda`` share one true eigenvector.

`T` must be entire inside the circle: the eigenvalue count is the winding
number of ``\det T``, which a pole would reduce.

**Errors**

- `"Maximum quadrature points reached"` if the quadrature has not converged
  after 65536 points. This usually means an eigenvalue lies on or very close
  to the unit circle.
- `SingularException` if `T(z)` is exactly singular at a quadrature point.
- An error if the count needs more than `Kmax` block rows,
  ``\lceil m / \min(n, \ell) \rceil > K_{\max}``. Increase `Kmax` or `nprobes`.
- An error if the winding number of ``\det T`` is negative (poles inside).

## `Beyn.countzeros`

```julia
Beyn.countzeros(T)
```

The number of eigenvalues of `T` inside the unit circle, with algebraic
multiplicity, from the winding number of ``\det T``. Assumes `T` is entire.

## Searching a different disc

Only the unit circle is built in. To search the disc ``|z - c| < \rho``,
rescale the problem and map the eigenvalues back:

```julia
c, ρ = 2.0, 0.5
sol = Beyn.solver(z -> T(c + ρ * z), 5)
λ = c .+ ρ .* sol.eigenvalues
```

The eigenvectors need no change.
