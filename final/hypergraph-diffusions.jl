##
include("../methods.jl")
include("../methods_hypergraph.jl")

H = _read_final_hypergraph("fauci-email-hypergraph.json")

##
A = project_hypergraph(H)
xy = igraph_layout(A; random=false)
# fill!(A.nzval, 1) # if you wish to make unweighted
G = (A = A, names = H.names, xy = xy,orgs = H.orgs)
## Compare Hypergraph Diffusions

# Include all the codes...
include("../include/LHQD-main/qdsfm-ppr.jl")
include("../include/LHQD-main/PageRank.jl")
include("../include/LHQD-main/common.jl")
include("../include/LHQD-main/local-hyper.jl")

##
function _make_diffusion_result(diffusion_type, x, names, topk)
  p = sortperm(x,rev=true)
  return string(diffusion_type) => (topk= names[p[1:topk]] .=> x[p[1:topk]], x, names)
end
results = Dict()
result_types = []
topk=15

##
kappa = 0.001  # value of kappa (sparsity regularization)
gamma = 0.1      # value of gamma (regularization on seed)
seedset = [nodeid(G,"conrad");]
deg = vec(sum(G.A,dims=2))
x = PageRank.acl_diffusion(G.A,deg,seedset,gamma, kappa)
push!(results, _make_diffusion_result("Graph-ACL", x, G.names, topk))
push!(result_types, "Graph-ACL" => "Sparse Seeded PageRank-Graph")

## Need to make hypergraph unweighted for this code...
H2 = unroll_weighted_hypergraph(H)
q = 2.0
rho = 0.0           # value of rho (KKT apprx-val)
L = LH.loss_type(q) # the loss-type, this is a 2-norm
GG = LH.graph(H2,1.0)
x_lh ,r,iter = LH.lh_diffusion(GG,seedset,gamma,kappa,rho,L,max_iters=100000)
push!(results, _make_diffusion_result("LHPR", x_lh, G.names, topk))
push!(result_types, "LHPR" => "Sparse Seeded PageRank-HyperGraph")

## QDSFM PPR
Ht = sparse(H2')
beta = 0.05
x_qd = QDSFMPageRank.qdsfmpr_ppr_euler(Ht, seedset, beta, 5000, 0.01)
push!(results, _make_diffusion_result("QDSFM", x_qd, G.names, topk))
push!(result_types, "QDSFM" => "Seeded PageRank-HyperGraph")

##
_write_score_table(results,result_types)
