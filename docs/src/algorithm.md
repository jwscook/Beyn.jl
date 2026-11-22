# Algorithm

## Contour moments

Pick ``\ell`` probe vectors, the columns of ``V \in \mathbb{C}^{n \times \ell}``,
and compute the contour moments of ``T^{-1}`` around the unit circle ``\Gamma``:

```math
A_p = \frac{1}{2\pi i} \oint_\Gamma z^p\, T(z)^{-1} V \,dz,
\qquad p = 0, 1, \dots, 2K_{\max} - 1 .
```

All ``2K_{\max}`` moments come from the same solves ``T(z)^{-1}V`` at the
quadrature points, so storing more of them costs memory and a few vector
updates per point, not more factorizations. The default is ``K_{\max} = 4``.

By Keldysh's theorem, if the ``k`` eigenvalues inside ``\Gamma`` are simple, then

```math
A_0 = X W^H V, \qquad A_1 = X \Lambda W^H V,
```

where the columns of ``X`` are the right eigenvectors, ``\Lambda`` is the
diagonal matrix of eigenvalues, and ``W`` holds the left eigenvectors.

## Counting the eigenvalues

The number ``m`` of eigenvalues inside ``\Gamma``, counted with algebraic
multiplicity, is found by the argument principle:

```math
m = \frac{1}{2\pi} \Delta_\Gamma \arg \det T(z),
```

the winding number of ``\det T`` around 0. This assumes ``T`` is **entire**
(analytic inside ``\Gamma``): a pole of ``\det T`` inside the contour would
be subtracted from the count. The phase of ``\det T(z)`` comes free from the
LU factorization at each quadrature point. The count sums the phase steps
between neighbouring points, so the solver bisects the angular grid until no
step exceeds ``\pi/4``. Each extra point costs one LU factorization.

## Block Hankel matrices

For ``K`` block rows, form the ``Kn \times K\ell`` block Hankel matrices

```math
H_0 = \begin{pmatrix}
A_0 & A_1 & \cdots & A_{K-1} \\
A_1 & A_2 & \cdots & A_K \\
\vdots & & & \vdots \\
A_{K-1} & A_K & \cdots & A_{2K-2}
\end{pmatrix}, \qquad
H_1 = \begin{pmatrix}
A_1 & \cdots & A_K \\
\vdots & & \vdots \\
A_K & \cdots & A_{2K-1}
\end{pmatrix}.
```

``H_0`` must have rank ``m``, which needs ``K \min(n, \ell) \ge m``. The
solver starts at ``K = \lceil m / \min(n,\ell) \rceil`` and increases ``K``,
up to ``K_{\max}``, until ``H_0`` has at least ``m`` singular values above the
cutoff **and is not saturated**. Saturated means every singular value is above
the cutoff: the rank is then capped by the size of ``H_0``, not by the number
of eigenvalues, so ``H_0`` may not have captured them all. A spare singular
value at the noise level confirms the rank has levelled off at ``m``. The
``K`` where this first happens is ``K_{\text{eig}}``. For example, two
``2\times2`` Jordan blocks with ``n = \ell = 2`` give a saturated ``4\times4``
``H_0`` at ``K = 2``, so the solver moves to ``K = 3``.
With ``K = 1`` this is the original method,
``H_0 = A_0`` and ``H_1 = A_1``. Larger ``K`` handles more eigenvalues than
probes, and eigenvalues with Jordan blocks.

## Reducing to a small linear problem

Take the singular value decomposition ``H_0 = V_0 \Sigma_0 W_0^H``, keeping
the first ``m`` singular values. Then

```math
B = V_0^H H_1 W_0 \Sigma_0^{-1}
```

is an ``m \times m`` matrix with the same eigenvalues as ``T`` inside
``\Gamma``. If ``B s = \lambda s``, the first ``n`` rows of ``V_0 s`` are an
eigenvector of the original problem.

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

## Checking the rank

The rank is the count ``m`` from the argument principle. The singular values
of ``H_0`` are checked against the cutoff

```math
\max\left(\texttt{svdrtol}\cdot\sigma_1,\; \texttt{svdatol},\; \texttt{svdrtol}\cdot K s\right),
```

where ``s`` approximates ``\frac{1}{2\pi}\oint \|T(z)^{-1}V\|\,|dz|``. The
solver warns if, at the final ``K``:

- ``\sigma_m`` is below the cutoff even at ``K = K_{\max}``: the subproblem is
  ill-conditioned;
- every singular value is still above the cutoff at ``K = K_{\max}``
  (saturated): the count cannot be confirmed, so try a larger ``K_{\max}`` or
  more probes;
- ``\sigma_{m+1}`` is above the cutoff: ``T`` may not be entire, or the
  quadrature is inaccurate.

## Requirements

- ``T`` must be entire inside ``\Gamma`` (no poles), so that the winding
  number of ``\det T`` counts its eigenvalues.
- ``\lceil m / \min(n, \ell) \rceil \le K_{\max}``, otherwise the solver
  stops with an error.
- ``\ell`` must be at least the largest geometric multiplicity of an
  eigenvalue inside ``\Gamma``.
- ``T`` must not be singular on ``\Gamma`` itself.

## Reference

W.-J. Beyn, *An integral method for solving nonlinear eigenvalue problems*,
Linear Algebra and its Applications 436 (2012) 3839–3863.
