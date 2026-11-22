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
  o = similar(probes, complex(float(eltype(probes))))
  function addpoint!(acc, t)
    z = cis(t)
    ldiv!(o, lu!(f(z)), probes)
    c = z / 2π
    cz = c * z
    acc.A₀ .= muladd.(c, o, acc.A₀)
    acc.A₁ .= muladd.(cz, o, acc.A₁)
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
  V, S, W = svd!(A₀)

  r = ngoodsingularvalues(S, max(svdatol, svdrtol * scale), svdrtol)

  V0 = @view V[:, 1:r]
  W0 = @view W[:, 1:r]
  V0tA₁ = mul!(similar(A₁, r, size(A₁, 2)), V0', A₁)
  B = mul!(similar(A₁, r, r), V0tA₁, W0)
  rdiv!(B, Diagonal(@view S[1:r]))

  λ, eigvecs = eigen!(B)
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
