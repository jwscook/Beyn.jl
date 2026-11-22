module Beyn
using LinearAlgebra

# Running (unweighted) trapezoid sums for the moments A[p+1] = Aₚ, p = 0, …, P-1,
# plus ∑‖f(z)⁻¹V‖ |dz|/2π, which sets a round-off floor when an integral is ≈ 0.
mutable struct ReductionStruct{T<:AbstractMatrix}
  A::Vector{T}
  scale::Float64
end
ReductionStruct(A::AbstractMatrix, P::Integer) = ReductionStruct([zero(A) for _ in 1:P], 0.0)
nmoments(acc::ReductionStruct) = length(acc.A)

function Base.copyto!(dst::ReductionStruct, src::ReductionStruct)
  foreach(copyto!, dst.A, src.A)
  dst.scale = src.scale
  return dst
end

# ‖αx - βy‖ without allocating a temporary
function normdiff(x, y, α, β)
  s = zero(real(eltype(x)))
  @inbounds for i in eachindex(x, y)
    s += abs2(muladd(α, x[i], -β * y[i]))
  end
  return sqrt(s)
end

# Trapezoid rule on [a, b] for a periodic integrand, doubling N until converged.
# `addpoint!(acc, t)` adds the integrand at t into acc. The sums in acc are kept
# unweighted (the integral with N points is h * sum) and scaled by h on return.
function adaptivetrapz!(addpoint!::P, acc::ReductionStruct, a, b; atol=1e-12, rtol=1e-10,
    N=16, maxN=65536) where {P}
  prev = ReductionStruct(acc.A[1], nmoments(acc))
  h = (b - a) / N
  for k in 0:N-1
    addpoint!(acc, muladd(k, h, a))
  end
  while 2N <= maxN
    copyto!(prev, acc)
    for k in 0:N-1
      addpoint!(acc, muladd(k + 0.5, h, a))
    end
    # I_N = h * prev, I_2N = (h/2) * acc, measured over all moments at once
    err = sqrt(sum(i -> normdiff(acc.A[i], prev.A[i], h / 2, h)^2, eachindex(acc.A)))
    nI = h / 2 * sqrt(sum(x -> norm(x)^2, acc.A))
    # On the unit circle |zᵖ| = 1, so every moment's integrand has norm ‖f⁻¹V‖/2π,
    # and ‖[A₀ … A_{P-1}]‖ integrand = √P × that.
    tol = max(atol, rtol * nI, 100eps() * h / 2 * sqrt(nmoments(acc)) * acc.scale)
    # The trapezoid rule converges exponentially for smooth periodic integrands,
    # so the error of I_2N is ≈ err²/‖I‖ rather than err.
    if err <= tol || err^2 <= tol * nI
      foreach(x -> lmul!(h / 2, x), acc.A)
      acc.scale *= h / 2
      return acc, err, 2N
    end
    N *= 2
    h /= 2
  end
  error("Maximum quadrature points reached")
end

# Phase of det f(z), from an LU factorization so it cannot overflow.
# logabsdet(F) returns (log|det|, det/|det|); only the unit phase factor is used.
detphase(F::LU) = angle(last(logabsdet(F)))

# Moments Aₚ = (1/2πi)∮ zᵖ f(z)⁻¹ V dz on the unit circle, p = 0, …, P-1.
# With z = e^{it}, dz = i z dt, so the integrand of Aₚ is zᵖ f(z)⁻¹ V z / 2π.
# f(z) must return a fresh matrix each call: it is factorized in place.
# Also returns the phase of det f(z) at every quadrature node, as (t, phase) pairs,
# which `countzeros` turns into the number of eigenvalues inside the circle.
function contourmoments(f::F, probes, P::Integer; atol=0.0, rtol=sqrt(eps())) where {F}
  P >= 2 || throw(ArgumentError("need at least 2 moments, got P = $P"))
  # Work buffer for f(z)⁻¹ V, allocated once and reused at every quadrature point.
  # Complex even if the probes are real, because f(z) is complex on the contour.
  o = similar(probes, complex(float(eltype(probes))))
  phases = Tuple{Float64,Float64}[]
  # Adds the integrand at angle t into the running sums. Nothing here allocates
  # except f(z) itself, the LU pivot vector and the growth of `phases`.
  function addpoint!(acc, t)
    z = cis(t)                    # point on the unit circle, z = e^{it}
    # o = f(z) \ probes, without temporaries:
    #   lu!(M)          factorizes M in place (overwrites the fresh matrix from f)
    #                   instead of copying it as lu(M) or M \ probes would;
    #   ldiv!(o, F, V)  solves F * o = V, writing the result into the buffer o.
    LUf = lu!(f(z))
    ldiv!(o, LUf, probes)
    push!(phases, (t, detphase(LUf)))  # the factorization gives det f(z) for free
    # c = (1/2πi) * dz/dt with dz/dt = i z: the weight for A₀. Each further
    # moment has one more factor of z. The fused in-place update
    # acc.A[p] .= muladd.(c, o, acc.A[p]) is a single loop writing straight
    # into acc.A[p], with no temporary c * o array.
    c = z / 2π
    for Ap in acc.A
      Ap .= muladd.(c, o, Ap)
      c *= z
    end
    # |zᵖ| = 1, so the integrand of every moment has norm ‖o‖/2π. Summed up,
    # this approximates ∮‖f⁻¹V‖|dz|/2π and sets the round-off floor.
    acc.scale += norm(o) / 2π
    return acc
  end
  acc, _, _ = adaptivetrapz!(addpoint!, ReductionStruct(o, P), 0, 2π; rtol, atol) # must be anti-clockwise
  return acc.A, acc.scale, phases
end

# Number of zeros of det f inside the unit circle, by the argument principle:
# the winding number of det f(e^{it}) around 0 as t goes from 0 to 2π. This
# counts eigenvalues with algebraic multiplicity, and is correct only when f is
# entire (analytic) inside the circle: a pole of det f would subtract from it.
#
# `phases` are (t, arg det f(e^{it})) pairs at any set of angles in [0, 2π). The
# winding number is the sum of phase steps between neighbouring angles, each
# taken in (-π, π]. That is only right if no step really exceeds π, so while the
# largest step is above `maxstep` the angular grid is refined by bisection,
# which costs one LU factorization of f per new point.
function countzeros(f::F, phases::Vector{Tuple{Float64,Float64}};
    maxstep=π / 4, maxN=65536) where {F}
  isempty(phases) && throw(ArgumentError("no phase samples given"))
  pts = sort(phases; by=first)
  while true
    N = length(pts)
    # next angle after pts[i], with the last one wrapping around to the first + 2π
    tnext(i) = i == N ? pts[1][1] + 2π : pts[i+1][1]
    step(i) = rem2pi(pts[mod1(i + 1, N)][2] - pts[i][2], RoundNearest)
    if maximum(i -> abs(step(i)), 1:N) <= maxstep
      return round(Int, sum(step, 1:N) / 2π)
    end
    2N <= maxN || error("det f(z) varies too fast on the contour to count its zeros " *
                        "with $maxN points; an eigenvalue may be on or near the unit circle")
    for i in 1:N
      t = (pts[i][1] + tnext(i)) / 2
      push!(pts, (t, detphase(lu!(f(cis(t))))))
    end
    sort!(pts; by=first)
  end
end
countzeros(f::F; kwargs...) where {F} =
  countzeros(f, [(t, detphase(lu!(f(cis(t))))) for t in range(0, 2π; length=17)[1:16]]; kwargs...)

# Number of singular values above the cutoff (S is sorted largest first).
ngoodsingularvalues(S, cutoff) = something(findlast(>(cutoff), S), 0)

# Block Hankel matrix [A_{i+j+shift}], i, j = 0, …, K-1, with n×ℓ blocks.
# H₀ (shift = 0) uses A₀ … A_{2K-2}, H₁ (shift = 1) uses A₁ … A_{2K-1}.
function blockhankel(A::AbstractVector{<:AbstractMatrix}, K::Integer, shift::Integer)
  n, ℓ = size(A[1])
  H = similar(A[1], K * n, K * ℓ)
  for j in 0:K-1, i in 0:K-1
    copyto!(view(H, i*n+1:(i+1)*n, j*ℓ+1:(j+1)*ℓ), A[i+j+shift+1])
  end
  return H
end

# Smallest K in Kmin:Kmax for which H₀ has at least m good singular values AND
# is not saturated, falling back to Kmax. Saturated means every singular value
# clears the cutoff: the rank is then limited by the size of H₀,
# min(K n, K ℓ), rather than by the number of eigenvalues, so H₀ may not have
# captured them all. A spare singular value at the noise level confirms the rank
# has plateaued at the true count. Returns the SVD of H₀ at that K, with the
# cutoff and the number r of singular values above it.
function chooseKeig(A, m, Kmin, Kmax, scale, svdatol, svdrtol)
  for K in Kmin:Kmax
    # svd! may overwrite H₀ instead of copying it, which is fine: it is not used again.
    V, S, W = svd!(blockhankel(A, K, 0))
    # ‖H₀‖ grows like K × the norm of one moment, so the round-off floor
    # svdrtol * scale grows by K as well.
    cutoff = max(svdatol, svdrtol * S[1], svdrtol * K * scale)
    r = ngoodsingularvalues(S, cutoff)
    if (r >= m && r < length(S)) || K == Kmax
      return (; V, S, W, Keig=K, cutoff, r)
    end
  end
  # Only reached if Kmin > Kmax; makes the return type exactly the NamedTuple.
  throw(ArgumentError("empty Keig range $Kmin:$Kmax"))
end

"""
    Beyn.solver(f, probes::AbstractMatrix; Kmax=4, kwargs...)
    Beyn.solver(f, nprobes::Integer; Kmax=4, kwargs...)

Eigenvalues of the entire matrix function `f` inside the unit circle, by Beyn's
method with block Hankel matrices of up to `Kmax` block rows.

The moments A₀ … A_{2Kmax-1} are computed in a single contour quadrature. The
number of eigenvalues `m` inside the circle (with algebraic multiplicity) is
found from the winding number of det f, and `Keig` is the smallest number of
block rows with `Keig * min(n, ℓ) ≥ m` for which the `Keig × Keig` block Hankel
matrix H₀ has at least `m` singular values above the cutoff and is not saturated
(not every singular value above the cutoff); `Keig ≤ Kmax`.
"""
function solver(f::F, probes::AbstractMatrix;
    Kmax::Integer=4,
    svdatol=0.0, svdrtol=sqrt(eps()),
    quadatol=0.0, quadrtol=sqrt(eps())) where {F}
  Kmax >= 1 || throw(ArgumentError("Kmax must be at least 1, got $Kmax"))
  n, ℓ = size(probes)

  # All 2Kmax moments come from the same f(z)⁻¹V solves; higher moments only
  # add cheap axpys per quadrature point.
  A, scale, phases = contourmoments(f, probes, 2Kmax; atol=quadatol, rtol=quadrtol)

  # m = number of eigenvalues inside, counted with algebraic multiplicity.
  m = countzeros(f, phases)
  m < 0 && error("det f has a negative winding number ($m) on the unit circle, " *
                 "so f has poles inside it; f must be entire")
  T = eltype(A[1])
  if m == 0
    return (eigenvalues=T[], eigenvectors=similar(A[1], n, 0),
            subproblemsingularvalues=real(T)[], count=0, Keig=0)
  end

  # H₀ is (Keig n) × (Keig ℓ), so its rank, which must reach m, is at most
  # Keig min(n, ℓ). Start from the smallest Keig allowing that.
  Kmin = cld(m, min(n, ℓ))
  Kmin <= Kmax || error("found $m eigenvalues inside the contour, which needs at least " *
                        "Keig = $Kmin block rows with $ℓ probes and n = $n, but Kmax = $Kmax; " *
                        "increase Kmax or the number of probes")

  (; V, S, W, Keig, cutoff, r) = chooseKeig(A, m, Kmin, Kmax, scale, svdatol, svdrtol)
  if r < m
    @warn "Block Hankel matrix has numerical rank below the eigenvalue count $m " *
          "even with Keig = Kmax = $Kmax; eigenvalues may be inaccurate. " *
          "Try more probes or a larger Kmax." σₘ = S[m] cutoff
  elseif r == length(S)
    @warn "Block Hankel singular values are still saturated (all $r above the cutoff) " *
          "at Keig = Kmax = $Kmax, so the eigenvalue count $m cannot be confirmed. " *
          "Try more probes or a larger Kmax." σₘᵢₙ = S[end] cutoff
  elseif r > m
    @warn "Block Hankel matrix has $r > $m singular values above the cutoff; " *
          "f may not be entire inside the contour, or the quadrature is inaccurate." σₘ₊₁ = S[m+1] cutoff
  end

  # Keep the first m singular vectors. @view refers to the columns of V and W
  # without copying them, as V[:, 1:m] would.
  V0 = @view V[:, 1:m]
  W0 = @view W[:, 1:m]
  # B = V0' * H₁ * W0 * inv(Diagonal(S[1:m])), built in steps into preallocated arrays:
  #   mul!(C, X, Y)  computes C = X * Y in place, writing into C.
  #                  V0' is a lazy adjoint, so no transposed copy is made.
  H₁ = blockhankel(A, Keig, 1)
  V0tH₁ = mul!(similar(H₁, m, size(H₁, 2)), V0', H₁)   # m × Keig ℓ
  B = mul!(similar(H₁, m, m), V0tH₁, W0)               # m × m
  # rdiv!(B, D) does B = B / D in place. With D diagonal this divides column j
  # of B by S[j], applying Σ0⁻¹ without forming 1 ./ S or an inverse matrix.
  rdiv!(B, Diagonal(@view S[1:m]))

  # B is a small linear eigenproblem with the same eigenvalues as f inside the
  # contour. eigen! may overwrite B, which is no longer needed.
  λ, eigvecs = eigen!(B)
  # Eigenvectors of B live in the m-dimensional basis V0. The first block row of
  # H₀ is built from A₀ … A_{Keig-1}, whose columns lie in the span of the
  # eigenvectors of f, so the first n rows of V0 * s give an eigenvector x with
  # f(λ) x ≈ 0. (For a Jordan block, only one of the repeated λ has a true
  # eigenvector; the others correspond to generalized eigenvectors.)
  X = (@view V0[1:n, :]) * eigvecs
  return (eigenvalues=λ, eigenvectors=X, subproblemsingularvalues=S, count=m, Keig=Keig)
end

function solver(f::F, nprobes::Integer; kwargs...) where {F}
  n = size(f(1.0 + 0im), 1)
  probes = randn(ComplexF64, n, nprobes) # it sees A \ probes
  foreach(normalize!, eachcol(probes))
  return solver(f, probes; kwargs...)
end

end # module Beyn
