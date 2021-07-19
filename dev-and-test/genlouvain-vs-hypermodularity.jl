## Test HyperMoularity.LambdaLouvain with implicit vs. explicit graphs...
# Conclusion -- the raw LambdaLouvain call on the temporal modularity
# matrix directly is equivalent to another implementation of Louvain clustering.

include("../methods.jl")
G = _build_email_tofrom_graph(data;keepfauci=false,mindegree=2)
##
A = G.A
##
using HyperModularity
cc = HyperModularity.VanillaModularity(Float64.(A),1.0)
##
include("../include/Optimal_LambdaCC.jl")
##
# k = vec(sum(G.A,dims=2))
# B = sparse(A - (2/sum(A))*(k*k'))
# ccs = HyperModularity.LambdaLouvain(Float64.(B),zeros(size(B,1)),0.0)

##
k = vec(sum(G.A,dims=2))
B = sparse(A - (1/sum(A))*(k*k'))
ccs = HyperModularity.LambdaLouvain(Float64.(B),zeros(size(B,1)),0.0)

##
m1 = compute_modularity(A, vec(sum(A;dims=2)), cc[:,end])
##
function compute_modularity(A::SparseMatrixCSC,c)
    obj = 0
    m = sum(nonzeros(A)/2)
    n = size(A,1)
    # index c...
    C = sparse(1:n, c, 1, n, maximum(c))
    Crv = rowvals(C)
    rv = rowvals(A)
    nz = nonzeros(A)
    for ci=1:size(C,2)
      Cset = Set(Crv[nzrange(C,ci)]) # index the cluster
      for nzci in nzrange(C,ci)
        i = Crv[nzci]
        for nzi in nzrange(A,i)
          if rv[nzi] in Cset
            obj += nz[nzi]
          end
        end
      end
    end
    return obj
end
m2 = compute_modularity(B, cc[:,end])/sum(A)
##
m3 = compute_modularity(B, ccs[:,end])/sum(A)
