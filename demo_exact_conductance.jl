include("methods.jl")

## Get graph
G = _build_email_tofrom_graph(data; keepcc=true,
    maxset=5, keepfauci=false,mindegree = 0, emailweight = false) |> igraph_layout
drawgraph(G)

##
include("include/exact_conductance_jump.jl")
##
G = (G..., A = spones!(G.A - Diagonal(G.A)))
## Exact conductance
S, cond = exact_conductance(G.A)

## Plots
drawgraph(G)
drawset!(G,S,markershape = :square)

## compare against ncut
# This takes a lot longer!
include("include/Optimal_LambdaCC.jl")
Clus, Lams, ncut_S = exact_normalized_cut(G.A)
##
ncutval = compute_min_norm_cut(G.A,Clus[end-1])
##
##
ncutupper = cond/sum(G.A)*2

##
drawgraph(G)
drawset!(G,ncut_S,markershape = :square)
## See differences in sets
setdiff(setdiff(1:n,ncut_S),S)
##
## End demo

## Compare the conductance cut in various graphs
G = _build_email_tofrom_graph(data;keepfauci=false) |> igraph_layout

S, cond = exact_conductance(G.A)
drawgraph(G)
drawset!(G,S,markershape = :square)

##
Clus, Lams, ncut_S = exact_normalized_cut(Float64.(G.A),2*cond/sum(G.A))
##
drawgraph(G)
drawset!(G,ncut_S,markershape = :square)

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

##
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"), mindegree=2) |> igraph_layout
spones!(G.A)
S, cond = exact_conductance(G.A)
##
drawgraph(G)
drawset!(G,S,markershape = :square)
##
Clus, Lams, ncut_S = exact_normalized_cut(Float64.(G.A),2*cond/sum(G.A))
##
drawgraph(G)
drawset!(G,ncut_S,markershape = :square)
## Try tofrom graph with weights...
G = _build_email_tofrom_graph(data;keepfauci=false,keepcc=true,emailweight=true,mindegree=2) |> igraph_layout
##
S, cond = exact_conductance(G.A - Diagonal(G.A))

## Plots -- these are boring and trivial...
drawgraph(G)
drawset!(G,S,markershape = :square)
