include("methods.jl")

## Get graph
kf = false
ms = 5
mindeg = 0
G = _build_email_tofrom_graph(data; keepcc = true,
    maxset=ms, keepfauci=kf,mindegree = mindeg, emailweight = false) |> igraph_layout
drawgraph(G,markersize = 5)

##
include("include/Optimal_LambdaCC.jl")

# Make this unweighted, if desired
A = Float64.(G.A)
fill!(A.nzval, 1)

## Exact normalized cut and modularity
# May take a while...
Clus, Lams, ncut_S = exact_normalized_cut(A)

# first clustering will be the modularity solutions
c = Clus[1]

## Display modularity clustering
drawgraph(G,alpha = 0.5)
for i = 1:maximum(c)
    drawset!(G, findall(x->x==i,c),markersize = 8)
end

# and optimal normalized cut solution
drawset!(G,ncut_S;markershape = :square,markersize = 8)
