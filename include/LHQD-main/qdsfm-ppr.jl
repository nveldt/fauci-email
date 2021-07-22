module QDSFMPageRank
using SparseArrays, LinearAlgebra
"""
assumes that each column (edge) has at least one entry in it
"""
function compute_subgradient!(H::SparseMatrixCSC, d::Vector, x::Vector, g::Vector)
  fill!(g, 0)
  rowval = rowvals(H)
  @inbounds for e = 1:size(H,2)
    # find the nodes with the extremal scores on the current hyperedge
    vmin = rowval[nzrange(H,e)[1]]
    vmax = vmin
    for nzi in nzrange(H, e)
      v = rowval[nzi]
      if x[v] < x[vmin]/d[vmin]
        vmin = v
      end
      if x[v] > x[vmax]/d[vmax]
        vmax = v
      end
    end
    g[vmax] += 2*(x[vmax]/d[vmax]-x[vmin]/d[vmin])
    g[vmin] += 2*(x[vmin]/d[vmin]-x[vmax]/d[vmax])
  end
end
""" Estimate a QDSFM-PageRank vector using an Euler integration strategy.
H - the hypergraph, where each column is a hyperedge.
seeds - the vector of initial seeds.
"""
function qdsfmpr_ppr_euler(H::SparseMatrixCSC, seeds::Vector, beta::Real, niter::Integer, h::Real)
  n = size(H,1) # number of vertices
  d = vec(sum(H; dims=2))
  x = zeros(n)
  g = zeros(n)
  seedvol = sum(@view d[seeds])
  for v in seeds
    x[v] += d[v]/seedvol
  end
  svec = copy(x) # copy the seed vector
  svec .*= beta
  t1 = time()
  for t=1:niter
    t2 = time()
    # @show t,t2-t1
    compute_subgradient!(H, d, x, g)
    g .*= -(1-beta) # compute -(1-beta)*subgrad of Laplacian
    g .+= svec # compute g = (beta*s) - (1-beta)*subgrad of Laplacian
    BLAS.axpy!(-beta, x, g) # compute g += (-beta)*x = beta*s - beta*x - (1-beta*subgrad of Laplacian
    g .*= h
    x .+= g
  end
  return x
end
end # end module
# include("common.jl")
# using Test
# @testset "QDSFM-PPR" begin
#   name = "TinyZoo"
#   M = matread("small-hypergraphs/"*name*".mat")
#   H = M["H"]
#   Ht = copy(H')
#   S = [1]
#   beta = 0.05
#   x = QDSFMPageRank.qdsfmpr_ppr_euler(Ht, S, beta, 2000, 0.01)
# end
