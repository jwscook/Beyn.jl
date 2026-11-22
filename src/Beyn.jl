module Beyn
using LinearAlgebra

# Running (unweighted) trapezoid sums for A₀ and A₁, plus ∑‖integrand‖,
# which sets a round-off floor when the integral is ≈ 0.
mutable struct ReductionStruct{T<:AbstractMatrix}
  A₀::T
  A₁::T
  scale::Float64
end
ReductionStruct(A::AbstractMatrix) = ReductionStruct(zero(A), zero(A), 0.0)

function Base.copyto!(dst::ReductionStruct, src::ReductionStruct)
  copyto!(dst.A₀, src.A₀)
  copyto!(dst.A₁, src.A₁)
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
  prev = ReductionStruct(acc.A₀)
  h = (b - a) / N
  for k in 0:N-1
    addpoint!(acc, muladd(k, h, a))
  end
  while 2N <= maxN
    copyto!(prev, acc)
    for k in 0:N-1
      addpoint!(acc, muladd(k + 0.5, h, a))
    end
    # I_N = h * prev, I_2N = (h/2) * acc
    err = hypot(normdiff(acc.A₀, prev.A₀, h / 2, h), normdiff(acc.A₁, prev.A₁, h / 2, h))
    nI = h / 2 * hypot(norm(acc.A₀), norm(acc.A₁))
    tol = max(atol, rtol * nI, 100eps() * h / 2 * acc.scale)
    # The trapezoid rule converges exponentially for smooth periodic integrands,
    # so the error of I_2N is ≈ err²/‖I‖ rather than err.
    if err <= tol || err^2 <= tol * nI
      lmul!(h / 2, acc.A₀)
      lmul!(h / 2, acc.A₁)
      acc.scale *= h / 2
      return acc, err, 2N
    end
    N *= 2
    h /= 2
  end
  error("Maximum quadrature points reached")
end

# A₀ = (1/2πi)∮ f(z)⁻¹ V dz, A₁ = (1/2πi)∮ z f(z)⁻¹ V dz on the unit circle.
# With z = e^{it}, dz = i z dt, so the integrand is f(z)⁻¹ V z / 2π.
# f(z) must return a fresh matrix each call: it is factorized in place.
function calculateA₀A₁(f::F, probes; atol=0.0, rtol=sqrt(eps())) where {F}
  # Work buffer for f(z)⁻¹ V, allocated once and reused at every quadrature point.
  # Complex even if the probes are real, because f(z) is complex on the contour.
  o = similar(probes, complex(float(eltype(probes))))
  # Adds the integrand at angle t into the running sums. Nothing here allocates
  # except f(z) itself and the LU pivot vector.
  function addpoint!(acc, t)
    z = cis(t)                    # point on the unit circle, z = e^{it}
    # o = f(z) \ probes, without temporaries:
    #   lu!(M)          factorizes M in place (overwrites the fresh matrix from f)
    #                   instead of copying it as lu(M) or M \ probes would;
    #   ldiv!(o, F, V)  solves F * o = V, writing the result into the buffer o.
    ldiv!(o, lu!(f(z)), probes)
    c = z / 2π                    # (1/2πi) * dz/dt with dz/dt = i z: the weight for A₀
    cz = c * z                    # A₁ has an extra factor of z
    # Fused in-place updates, acc.A₀ += c * o and acc.A₁ += cz * o:
    #   the dots make each line a single loop over the elements that writes
    #   straight into acc.A₀ / acc.A₁, so no temporary c * o array is made;
    #   muladd(x, y, w) computes x * y + w, as one fused multiply-add where
    #   the hardware supports it.
    acc.A₀ .= muladd.(c, o, acc.A₀)
    acc.A₁ .= muladd.(cz, o, acc.A₁)
    # Norm of this point's contribution to [A₀ A₁]: ‖[c o, cz o]‖ = |c| √(1 + |z|²) ‖o‖.
    # Summed up, this approximates ∫‖integrand‖ and sets the round-off floor.
    acc.scale += abs(c) * sqrt(1 + abs2(z)) * norm(o)
    return acc
  end
  acc, _, _ = adaptivetrapz!(addpoint!, ReductionStruct(o), 0, 2π; rtol, atol) # must be anti-clockwise
  return acc.A₀, acc.A₁, acc.scale
end

function ngoodsingularvalues(S, svdatol, svdrtol)
  return something(findlast(s -> s > max(S[1] * svdrtol, svdatol), S), 0)
end
function solver(f::F, probes;
    svdatol=0.0, svdrtol=sqrt(eps()),
    quadatol=0.0, quadrtol=sqrt(eps())) where {F}

  A₀, A₁, scale = calculateA₀A₁(f, probes, atol=quadatol, rtol=quadrtol)
  # A₀ = V * Diagonal(S) * W'. svd! may overwrite A₀ instead of copying it,
  # which is fine because A₀ is not used again.
  V, S, W = svd!(A₀)

  # Rank r = number of eigenvalues found. The svdrtol * scale term is an absolute
  # floor so that a pure round-off A₀ (no eigenvalues inside) gives r = 0.
  r = ngoodsingularvalues(S, max(svdatol, svdrtol * scale), svdrtol)

  # Keep the first r singular vectors. @view refers to the columns of V and W
  # without copying them, as V[:, 1:r] would.
  V0 = @view V[:, 1:r]
  W0 = @view W[:, 1:r]
  # B = V0' * A₁ * W0 * inv(Diagonal(S[1:r])), built in steps into preallocated arrays:
  #   mul!(C, X, Y)  computes C = X * Y in place, writing into C.
  #                  V0' is a lazy adjoint, so no transposed copy is made;
  #   similar(A₁, m, k) allocates an uninitialized m×k array of A₁'s element type.
  V0tA₁ = mul!(similar(A₁, r, size(A₁, 2)), V0', A₁)   # r×ℓ = V0' * A₁
  B = mul!(similar(A₁, r, r), V0tA₁, W0)               # r×r = (V0' * A₁) * W0
  # rdiv!(B, D) does B = B / D in place. With D diagonal this divides column j
  # of B by S[j], applying Σ0⁻¹ without forming 1 ./ S or an inverse matrix.
  rdiv!(B, Diagonal(@view S[1:r]))

  # B is a small linear eigenproblem with the same eigenvalues as f inside the
  # contour. eigen! may overwrite B, which is no longer needed.
  λ, eigvecs = eigen!(B)
  # Eigenvectors of B live in the r-dimensional basis V0; V0 * s maps them back
  # to eigenvectors x of the original problem, with f(λ) x ≈ 0.
  return (eigenvalues=λ, eigenvectors=V0 * eigvecs, subproblemsingularvalues=S)
end

function solver(f::F, nprobes::Integer;
    svdatol=0.0, svdrtol=sqrt(eps()),
    quadatol=0.0, quadrtol=sqrt(eps())) where {F}
  n = size(f(1.0 + 0im), 1)
  probes = randn(ComplexF64, n, nprobes) # it sees A \ probes
  foreach(normalize!, eachcol(probes))
  return solver(f, probes; svdatol=svdatol, svdrtol=svdrtol, quadatol=quadatol, quadrtol=quadrtol)
end

end # module Beyn
