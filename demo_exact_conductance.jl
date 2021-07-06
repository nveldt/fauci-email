include("methods.jl")

## Get graph
kf = false
ms = 5
mindeg = 2
G = _build_email_tofrom_graph(data; keepcc = true,
    maxset=ms, keepfauci=kf,mindegree = mindeg, emailweight = false) |> igraph_layout
drawgraph(G)

##
include("include/exact_conductance_jump.jl")

# Make this unweighted, if desired
n = size(A,1)
A = Float64.(G.A)
for i = 1:n
    A[i,i] = 0.0
end
dropzeros!(A)
fill!(A.nzval, 1)

## Exact conductance
S, cond = exact_conductance(A)

## Plots
drawgraph(G)
drawset!(G,S,markershape = :square)

## compare against ncut
# This takes a lot longer!
Clus, Lams, ncut_S = exact_normalized_cut(A,2*cond)


## Compare the conductance cut in various graphs
G = _build_email_tofrom_graph(data;keepfauci=false) |> igraph_layout
S, cond = exact_conductance(G.A)
drawgraph(G)
drawset!(G,S,markershape = :square)
## Is this the min collins,conrad st cut?
S_collins_conrad = stcut(G, "collins", "conrad")
setdiff(S_collins_conrad, setdiff(1:size(G.A,1),S)) # well, we have a subset
## this one is interesting...
# It's super large! But it is really half-vol. Robert Kadlac makes some
# very dense projections in the hypergraph as the emails involve
# a lot of people are are cited frequently.
G = _build_email_hypergraph_projection(data;keepfauci=false) |> igraph_layout
S, cond = exact_conductance(G.A)
drawgraph(G)
drawset!(G,S,markershape = :square)
## Remove weights
G = _build_email_hypergraph_projection(data;keepfauci=false) |> igraph_layout
spones!(G.A)
S, cond = exact_conductance(G.A)
drawgraph(G)
drawset!(G,S,markershape = :square)

## Let's try it with weights that are designed to
G = _build_email_hypergraph_projection(data;keepfauci=false,emailweight=true,mindegree=2) |> igraph_layout
S, cond = exact_conductance(G.A)
drawgraph(G)
drawset!(G,S,markershape = :square)
