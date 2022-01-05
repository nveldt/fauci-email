## Load the big graph
using MAT
data = matread("hypergraphs/VegasRestaurants_H.mat")
VegasH = data["H"]

## Try ACL with an implicit clique expansion.
##
module ImplicitClique
import SparseArrays: sparse
using SparseArrays, LinearAlgebra
struct ImplicitCliqueExpansion{T <: Real, Ti <: Integer}
  H::SparseMatrixCSC{T, Ti}
  d::Vector{T}
  c::Vector{T} # correction term
end
function _sum_of_eachrows_nonzeros_squared(A::SparseMatrixCSC)
  v = zeros(size(A,2))
  #rowval = rowvals(A)
  nzval = nonzeros(A)
  for j=1:size(A,2)
    for nzi in nzrange(A,j)
      v[j] += nzval[nzi]^2
    end
  end
  return v
end
"""
- `H` The hypergraph incidence matrix where each row is a hyper-edge.
"""
function ImplicitCliqueExpansion(H::SparseMatrixCSC)
  c = _sum_of_eachrows_nonzeros_squared(H)
  d = H'*(H*ones(size(H,2))) .- c # node degrees
  T = eltype(H)
  Ti = eltype(rowvals(H))
  return ImplicitCliqueExpansion{T,Ti}(H, d, c)
end
function ImplicitWeightedCliqueExpansion(H::SparseMatrixCSC)
  ImplicitCliqueExpansion(Diagonal(1.0./sqrt.(vec(sum(H;dims=2))))*H)
end
SparseArrays.sparse(G::ImplicitCliqueExpansion) = G.H'*G.H - Diagonal(G.c)

function _inverse_scale_vector!(r::Vector, d::Vector)
  for i = 1:length(r)
    r[i] /= d[i]
  end
end
function orig_acliter(G::ImplicitCliqueExpansion,
    alpha::Real,
    seed,
    kappa::Real)
  n = size(G.H,2)
  x = zeros(n)
  r = zeros(n)
  d = G.d
  svol = 0.0
  for s in seed
      svol += d[s]
  end
  for s in seed
      r[s] = d[s]/svol
  end
  tau = kappa/(alpha*sum(@view d[seed]))
   # should use more efficient data structure
  iter = 0
  while true
    # if iter % 1 == 0
    #   @show iter,sum(r)
    # end
    iter += 1
    T = findall(r .> tau*d)
    if length(T) == 0
      break
    end
    rT = zeros(n)
    rT[T] = r[T]-0.5*tau*d[T]
    r .-= rT
    x .+= (1-alpha)*rT
    #rT ./= d
    _inverse_scale_vector!(rT, d)
    r .+= alpha*(G.H'*(G.H*(rT)))
    r .-= alpha*G.c.*rT
    # r .+= alpha*A[:,T]*(rt./d[T]) # NOT DONE
    # r .-= rT
  end
  return x, r
end

function acliter(G::ImplicitCliqueExpansion,
    alpha::Real,
    seed,
    kappa::Real)

  n = size(G.H,2)
  x = zeros(n)
  r = zeros(n)
  z = zeros(size(G.H,1))
  d = G.d
  H = G.H
  svol = 0.0
  for s in seed
      svol += d[s]
  end
  for s in seed
      r[s] = d[s]/svol
  end
  tau = kappa/(alpha*sum(@view d[seed]))

  rowval = rowvals(G.H)
  nzval = nonzeros(G.H)
   # should use more efficient data structure
  iter = 0
  while true
    iter += 1

    nexcess = 0
    fill!(z, 0) # reset z

    @inbounds for v=1:n
      # T = findall(r .> tau*d)
      if r[v] > tau*d[v]
        #rT[T] = r[T]-0.5*tau*d[T]
        #r .-= rT
        nexcess += 1
        rv = r[v] - 0.5*tau*d[v]
        x[v] += (1-alpha)*(rv)
        r[v] -= rv
        # everything with scaled resid... TODO, precompute 1./d
        rv = rv./d[v]
        r[v] -= alpha*G.c[v]*rv # extra correction...  #r .-= alpha*G.c.*rT
        # compute alpha*G.H*rT
        for nzi = nzrange(H, v)
          u = rowval[nzi]
          z[u] += alpha*nzval[nzi]*rv
        end
      end
    end
    if nexcess == 0
      break
    end

    # r .+= G.H'*z
    @inbounds for v=1:n
      for nzi = nzrange(H, v)
        r[v] += nzval[nzi]*z[rowval[nzi]]
      end
    end
  end
  return x, r
end

function _acliter_simple(A::AbstractMatrix, alpha::Real, seed, kappa::Real, d::Vector=vec(sum(A;dims=1)))
  n = size(A,1)
  x = zeros(n)
  r = zeros(n)
  svol = 0.0
  for s in seed
      svol += d[s]
  end
  for s in seed
      r[s] = d[s]/svol
  end
  tau = kappa/(alpha*sum(@view d[seed]))
   # should use more efficient data structure
  while true
    T = findall(r .> tau*d)
    if length(T) == 0
      break
    end
    rt = r[T] - 0.5*tau*d[T]        # estimate excess residual
    x[T] = x[T] .+ (1-alpha)*rt # add excess residual to solution
    r .+= alpha*A[:,T]*(rt./d[T]) # NOT DONE
    r[T] .-= rt
  end
  return x, r
end
end # end module
using Test
@testset "ImplicitACL" begin
  #=using NearestNeighbors
  function random_graph_model(n::Integer, m::Integer)
    xy = randn(2,n)
  end=#
  using MAT, SparseArrays, LinearAlgebra
  data = matread("small-hypergraphs/TinyZoo.mat")
  H = data["H"]
  G = ImplicitClique.ImplicitCliqueExpansion(H)
  A = sparse(G)
  d = vec(sum(A;dims=2))
  @test d ≈ G.d
  alpha = 0.85
  seed = [1]
  tau = 1e-4
  x,xr = ImplicitClique._acliter_simple(A, alpha, seed, tau)
  y,yr = ImplicitClique.orig_acliter(G, alpha, seed, tau)
  z,zr = ImplicitClique.acliter(G, alpha, seed, tau)
  @test x ≈ y
  @test x ≈ z
  @test y ≈ z
  @test all(xr .>=0)
  @test all(yr .>=0)
  @test all(zr .>=0)
  @test all(xr .<= tau*d)
  @test all(yr .<= tau*d)
  @test all(zr .<= tau*d)
  @test norm(y-z) <= eps(1.0)
  println("")
  println("Norms")
  println(norm(y-z))
  println(norm(x-y))
  #@code_warntype(ImplicitClique.fast_acliter(G, alpha, seed, tau))
end
##
#ImplicitClique._sum_of_eachrows_nonzeros_squared(data["H"])
using BenchmarkTools
function perftest()

  G1 = ImplicitClique.ImplicitCliqueExpansion(VegasH)
  G2 = ImplicitClique.ImplicitWeightedCliqueExpansion(VegasH)
  alpha = 0.9
  seed = 1
  tau = 1e-4
  println("")
  println("Original Algs")
  @btime y1,y2r = ImplicitClique.orig_acliter($G1, $alpha, $seed, $tau) seconds=1
  @btime y2,y2r = ImplicitClique.orig_acliter($G2, $alpha, $seed, $tau) seconds=1
  println("Fast Algs")
  @btime fy1,fy1r = ImplicitClique.acliter($G1, $alpha, $seed, $tau) seconds=1
  @btime fy2,fy2r = ImplicitClique.acliter($G2, $alpha, $seed, $tau) seconds=1

end
perftest()
