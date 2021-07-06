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
Clus, Lams, ncut_S = exact_normalized_cut(A,2*cond)

