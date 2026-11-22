# API

## `Beyn.solver`

```julia
Beyn.solver(T, nprobes::Integer; contour=Beyn.Ellipse(), kwargs...)
Beyn.solver(T, probes::AbstractMatrix; contour=Beyn.Ellipse(), kwargs...)
```

Find the eigenvalues of ``T(\lambda)x = 0`` inside `contour`.

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
| `contour`  | `Ellipse()`   | an `Ellipse` or `Rectangle` to search (see [Contours](#Contours)); default the unit circle |
| `Kmax`     | `4`           | maximum block rows of the Hankel matrices; moments ``A_0 \dots A_{2K_{\max}-1}`` are stored |
| `quadrtol` | `sqrt(eps())` | relative tolerance of the contour quadrature             |
| `quadatol` | `0.0`         | absolute tolerance of the contour quadrature             |
| `svdrtol`  | `sqrt(eps())` | relative cutoff for the singular values of ``H_0``       |
| `svdatol`  | `0.0`         | absolute cutoff for the singular values of ``H_0``       |

**Returns** a named tuple with:

| field                      | contents                                                       |
|----------------------------|----------------------------------------------------------------|
| `eigenvalues`              | vector of eigenvalues inside the contour                       |
| `eigenvectors`             | matrix whose `i`th column is the eigenvector for `eigenvalues[i]` |
| `subproblemsingularvalues` | all singular values of ``H_0``, largest first                  |
| `count`                    | eigenvalues inside, with multiplicity, from the argument principle |
| `Keig`                     | block rows used, ``\le K_{\max}`` (0 if `count == 0`)           |
| `contour`                  | the `Ellipse` that was searched                                |

The eigenvectors are not normalized. For an eigenvalue with a Jordan block,
the repeated copies of ``\lambda`` share one true eigenvector.

`T` must be entire inside the contour: the eigenvalue count is the winding
number of ``\det T``, which a pole would reduce.

**Errors**

- `"Maximum quadrature points reached"` if the quadrature has not converged
  after 65536 points. This usually means an eigenvalue lies on or very close
  to the contour.
- `SingularException` if `T(z)` is exactly singular at a quadrature point.
- An error if the count needs more than `Kmax` block rows,
  ``\lceil m / \min(n, \ell) \rceil > K_{\max}``. Increase `Kmax` or `nprobes`.
- An error if the winding number of ``\det T`` is negative (poles inside).

## `Beyn.countzeros`

```julia
Beyn.countzeros(T; contour=Beyn.Ellipse())
```

The number of eigenvalues of `T` inside `contour`, with algebraic
multiplicity, from the winding number of ``\det T``. Assumes `T` is entire.

## Contours

### `Beyn.Ellipse`

```julia
Beyn.Ellipse(centre=0, a=1, b=a; angle=0)
```

The ellipse ``z(t) = c + e^{i\theta}(a\cos t + i b \sin t)``, traversed
anti-clockwise, with centre ``c``, semi-axis ``a`` along the direction
``\theta`` = `angle` and semi-axis ``b`` perpendicular to it. `Ellipse()` is
the unit circle and `Ellipse(c, r)` the circle ``|z - c| = r``. The semi-axes
must be positive.

### `Beyn.Rectangle`

```julia
Beyn.Rectangle(z1, z2)
Beyn.Ellipse(rectangle::Beyn.Rectangle; margin=0.0)
```

An axis-aligned rectangle with opposite corners `z1` and `z2`; it must have
positive width and height. When passed as `contour`, it is replaced by
`Beyn.Ellipse(rectangle)`: the axis-aligned, least-area ellipse through its
corners, with semi-axes ``\sqrt{2}`` times the half-width and half-height. The
corners lie on that ellipse, so to keep the contour away from an eigenvalue near
a corner pass `contour=Beyn.Ellipse(rectangle; margin=0.1)`, which enlarges both
semi-axes by 10%.

Every eigenvalue inside the ellipse is returned, including those outside the
rectangle. To keep only the ones in the rectangle:

```julia
rect = Beyn.Rectangle(0.5 - 0.1im, 1.2 + 0.2im)
sol = Beyn.solver(T, 3; contour=rect)
keep = in(rect).(sol.eigenvalues)
λ, X = sol.eigenvalues[keep], sol.eigenvectors[:, keep]
```

`in(z, ellipse)` and `in(z, rectangle)` test whether `z` is strictly inside.
