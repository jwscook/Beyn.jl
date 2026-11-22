# Algorithm

## Contour moments

Pick ``\ell`` probe vectors, the columns of ``V \in \mathbb{C}^{n \times \ell}``,
and compute the first two contour moments of ``T^{-1}`` around the unit circle
``\Gamma``:

```math
A_0 = \frac{1}{2\pi i} \oint_\Gamma T(z)^{-1} V \,dz, \qquad
A_1 = \frac{1}{2\pi i} \oint_\Gamma z\, T(z)^{-1} V \,dz .
```

By Keldysh's theorem, if the ``k`` eigenvalues inside ``\Gamma`` are simple, then

```math
A_0 = X W^H V, \qquad A_1 = X \Lambda W^H V,
```

where the columns of ``X`` are the right eigenvectors, ``\Lambda`` is the
diagonal matrix of eigenvalues, and ``W`` holds the left eigenvectors.

## Reducing to a small linear problem

Take the singular value decomposition ``A_0 = V_0 \Sigma_0 W_0^H``, keeping
only the ``k`` singular values above the cutoff. Then

```math
B = V_0^H A_1 W_0 \Sigma_0^{-1}
```

is a ``k \times k`` matrix with the same eigenvalues ``\Lambda``. If
``B s = \lambda s``, then ``x = V_0 s`` is an eigenvector of the original
problem.

## Quadrature

With ``z = e^{it}``, ``dz = i z\,dt``, so

```math
A_0 = \int_0^{2\pi} T(e^{it})^{-1} V \,\frac{e^{it}}{2\pi}\,dt .
```

The integrand is smooth and periodic, so the trapezoid rule converges
exponentially. The solver starts with 16 points and doubles the number of
points, reusing all previous evaluations, until the estimated error meets the
tolerance. Because convergence is exponential, the error after a doubling is
roughly the square of the difference between the last two estimates, and the
solver uses that to stop one doubling earlier. The limit is 65536 points.

Each quadrature point costs one evaluation of ``T`` and one LU factorization.

## Choosing the rank

The rank ``k`` is the number of singular values of ``A_0`` above

```math
\max\left(\texttt{svdrtol}\cdot\sigma_1,\; \texttt{svdatol},\; \texttt{svdrtol}\cdot s\right),
```

where ``s`` approximates ``\oint \|T(z)^{-1}V\|\,|dz|``. The last term
matters when no eigenvalues are inside the contour: ``A_0`` is then pure
round-off, and a cutoff relative to ``\sigma_1`` alone would keep that noise.

## Requirements

- The number of eigenvalues inside ``\Gamma`` must satisfy ``k \le \ell \le n``.
- Eigenvalues must be simple. Jordan blocks need higher moments
  ``A_p = \frac{1}{2\pi i}\oint z^p T(z)^{-1} V\,dz``, which are not implemented.
- ``T`` must not be singular on ``\Gamma`` itself.

## Reference

W.-J. Beyn, *An integral method for solving nonlinear eigenvalue problems*,
Linear Algebra and its Applications 436 (2012) 3839–3863.
