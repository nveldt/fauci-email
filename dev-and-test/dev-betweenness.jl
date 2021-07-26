## dev betweeness centrality
include("../methods.jl")


##
G = _read_final_with_products("fauci-email-tofrom-5.json")
G = (G..., groups=G.products.simple.modularity)

##
function igraph_betweenness(A::SparseMatrixCSC{T}) where T
  ei,ej,ew = findnz(A)
  edgelist = [(ei[i]-1,ej[i]-1) for i = 1:length(ei)]

  nverts = size(A)
  G = igraph.Graph(nverts, edges=edgelist, directed=false)
  bc = G.betweenness()
  return bc
end
bc = igraph_betweenness(G.A)

## Test this with a drawing
G = (G..., groups=G.products.simple.modularity)
drawgraph(G, linecolor=:black, size=(550,550))
scatter!(G.xy[:,1],G.xy[:,2], markersize=log.(bc.+1).+2, color=G.groups,
  markerstrokewidth=1)
p = sortperm(bc;rev=true)
for i=p[1:15]
  #annotate!(G.xy[i,1],G.xy[i,2], G.names[i])
  showlabel!(G, G.names[i], 7, :left, theme_palette(:default)[G.groups[i]]; offset=3, fontargs=(;rotation=-rand(-25:25)))
end
plot!()
