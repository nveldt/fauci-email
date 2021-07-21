#=
This was just me figuring out how to write igraph_layout so we could get
non-random layouts.

These are initailied with circle layouts...
=#

include("methods.jl")
##
G = _build_email_repliedto_graph(data; keepfauci=false)
## Figure out how to get non-random layouts from igraph
function igraph_graph(A::SparseMatrixCSC{T}) where T
    ei,ej,ew = findnz(A)
    edgelist = [(ei[i]-1,ej[i]-1) for i = 1:length(ei)]
    nverts = size(A)
    G = igraph.Graph(nverts, edges=edgelist, directed=false)
    return G
end
pG = igraph_graph(G.A)
##
C = pG.layout_circle()
## Using the same initialization doesn't work??
xy1 = pG.layout_fruchterman_reingold(seed=collect(C))
##
xy2 = pG.layout_fruchterman_reingold(seed=collect(C))
##
xy2[1], xy1[1]
## Try setting the pyrandom generator
py"""import random
random.seed(0)"""
xy1 = pG.layout_fruchterman_reingold()
##
py"""import random
random.seed(0)"""
xy2 = pG.layout_fruchterman_reingold()
##
xy2[1], xy1[1]
##
include("../methods.jl")
xy1 = igraph_layout(G; random=false)
xy2 = igraph_layout(G; random=false)
norm(xy1.xy .- xy2.xy)
