include("methods_hypergraph.jl")
include("methods.jl")

## Build hypergraph
kf = false
ms = 10
mindeg = 5
parts = ("sender","recipients")
H = _build_email_hypergraph(data;hyperedgeparts=parts,maxset=ms, keepfauci=kf,mindegree = mindeg)
s = nodeid(H,"collins");
t = nodeid(H,"conrad");
S_hyp = hyperstcut(H.H, s,t,H.weights);
m,n = size(H.H)
if in(s,S_hyp)
    S_hyp = setdiff(1:n,S_hyp)
end


## Hypergraph projection into graph
# (careful, this is different in subtle ways from calling
# _build_email_hypergraph_projection, because of weighting and degree stuff)

A = project_hypergraph(H)
xy = igraph_layout(A)
# fill!(A.nzval, 1) # if you wish to make unweighted
G = (A = A, names = H.names, xy = xy,orgs = H.orgs)
# S = stcut(G, "collins", "conrad")
S = stcut(A, s,t; smallside= false)
if in(s,S)
    S = setdiff(1:n,S)
end

## Plot different node sets
JustHype = setdiff(S_hyp,S)
JustGraph = setdiff(S,S_hyp)
both = intersect(S,S_hyp)
drawgraph(G)
drawset!(G, both,markersize = 5, color =:purple, label = "both")
drawset!(G, JustHype,markersize = 5, color =:blue, label = "just-hyper")
drawset!(G, JustGraph,markersize = 5, color =:red, label = "just-graph")
drawset!(G,[s],markersize = 10)
drawset!(G,[t],markersize = 10)

savefig("bothstcut_mindeg$(mindeg)_maxset$(ms)_fauci_$(kf).pdf")


## What are these nodes where sets differ?
special = [JustHype; JustGraph]
for j in special
    println(rpad("$(G.names[j])", 25))
end

## PageRank
include("include/LHQD-main/qdsfm-ppr.jl")
include("include/LHQD-main/PageRank.jl")
include("include/LHQD-main/common.jl")
include("include/LHQD-main/local-hyper.jl")

##
kappa = 0.001  # value of kappa (sparsity regularization)
gamma = 0.1      # value of gamma (regularization on seed)
seedset = [nodeid(G,"conrad");]
deg = vec(sum(G.A,dims=2))
x = PageRank.acl_diffusion(G.A,deg,seedset,gamma, kappa)

## Need to make hypergraph unweighted for this code...
H2 = unroll_weighted_hypergraph(H)

## LHPR
q = 2.0
rho = 0.0           # value of rho (KKT apprx-val)
L = LH.loss_type(q) # the loss-type, this is a 2-norm
GG = LH.graph(H2,1.0)
x_lh ,r,iter = LH.lh_diffusion(GG,seedset,gamma,kappa,rho,L,max_iters=10000)

## QDSFM PPR
Ht = sparse(H2')
beta = 0.05
x_qd = QDSFMPageRank.qdsfmpr_ppr_euler(Ht, S, beta, 5000, 0.01)

## Rankings
p = sortperm(x,rev = true)
p2 = sortperm(x_qd,rev = true)
p3 = sortperm(x_lh,rev = true)
thenames = G.names
top = 10
println("Graph PPR  \t LHPR \t QDSFM-PPR")
println("------------------------------------------")
for k = 1:top
    print(rpad("$(thenames[p[k]])", 25))
    print(rpad("$(thenames[p3[k]])", 25))
    println(rpad("$(thenames[p2[k]])", 25))
end
println("------------------------------------------")
for k = 1:top
    println("$(x[p[k]])  \t $(x_lh[p3[k]]) \t $(x_qd[p2[k]])")
end


## Plot ppr top 10
top = 10
lh_clus = p3[1:top]
qd_clus = p2[1:top]
graph_clus = p[1:top]

drawgraph(G)
drawset!(G, graph_clus,markersize = 15,markercolor = :red)
drawset!(G, lh_clus,markersize = 10,color = :blue)
drawset!(G, qd_clus,markersize = 5, color = :green)

savefig("top$(top)_PageRank_mindeg$(mindeg)_maxset$(ms)_fauci_$(kf).pdf")
