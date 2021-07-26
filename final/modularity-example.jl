##
include("../methods.jl")
##
G = _read_final_with_products("fauci-email-tofrom-5.json")
G = (G..., groups=G.products.simple.modularity, xy=G.products.simple.xy)

drawgraph(G, size=(350,350))
bc = igraph_betweenness(G.A)
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1, hover=G.names[dp])
##
showlabel!(G,"collins",7)
showlabel!(G,"conrad",7)
showlabel!(G,"giroir",7)
showlabel!(G,"redfield",7)
showlabel!(G,"farrar",7)
showlabel!(G,"niaid",7)
showlabel!(G,"routh",7)
showlabel!(G,"niaid news",7)
showlabel!(G,"kadlec",7)
showlabel!(G,"new england",7)
showlabel!(G,"kaplan",7)
showlabel!(G,"giroir",7)
showlabel!(G,"birx",7)
##
showlabel!(G,"awwad",7)
showlabel!(G,"abutaleb",7)
showlabel!(G,"beigel",7)
showlabel!(G,"cabezas",7)
showlabel!(G,"auchin",7)

## Test this with a drawing
function _short_name(s::AbstractString)
  if occursin(",",s)
    return string(split(s,",")[1])
  elseif occursin(" ", s)
    return join(split(s, " "), "\n")
  end
end
drawgraph(G, linecolor=:black, size=(550,550))
# show groups ...
using LazySets
for gid = 1:maximum(G.groups)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.75) for i in findall(g->g==gid,G.groups)])
  plot!(h,alpha=0.1, color=gid, linewidth=0)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.5) for i in findall(g->g==gid,G.groups)])
  plot!(h,alpha=0.1, color=gid, linewidth=0)
end
# permute by bc size...
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1)
for i=dp[1:14]
  #annotate!(G.xy[i,1],G.xy[i,2], G.names[i])
  showlabel!(G,G.names[i], 9, :center,
      _darker(theme_palette(:default)[G.groups[i]], 0.66),
    #theme_palette(:default)[G.groups[i]]
    ; offset=0, textfunc=_short_name,
    #fontargs=(;rotation=-rand(-25:25))
    )
end
_darker(x::RGB, f::Real) = RGB(f*x.r, f*x.g, f*x.b)
# show other labels
for n in ["awwad", "abutaleb", "beigel", "cabezas", "niaid news"]
  nid = nodeid(G, n)
  showlabel!(G,n ,7, :center,
    _darker(theme_palette(:default)[G.groups[nid]], 0.5),
    textfunc=_short_name)
  end
plot!()
##
savefig("figures/modularity-tofrom.pdf")
## Try adding in CCs
G = _read_final_with_products("fauci-email-tofrom-cc-5.json")
G = (G..., groups=G.products.weighted.modularity, xy=G.products.weighted.xy)
drawgraph(G, size=(350,350))
bc = igraph_betweenness(G.A)
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1, hover=G.names[dp])

## Test this with a drawing
drawgraph(G, linecolor=:black, size=(550,550))
# show groups ...
using LazySets
for gid = 1:maximum(G.groups)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.75) for i in findall(g->g==gid,G.groups)])
  if gid > 16
    plot!(h,alpha=0.1, color=:black, linewidth=0)
  else
    plot!(h,alpha=0.1, color=gid, linewidth=0)
  end
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.5) for i in findall(g->g==gid,G.groups)])
  if gid > 16
    plot!(h,alpha=0.1, color=:black, linewidth=0)
  else
    plot!(h,alpha=0.1, color=gid, linewidth=0)
  end
end
# permute by bc size...
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1)
for i=dp[1:16]
  #annotate!(G.xy[i,1],G.xy[i,2], G.names[i])
  tcolor = G.groups[i] > 16 ? :black : _darker(theme_palette(:default)[G.groups[i]], 0.5)
  showlabel!(G,G.names[i], 9, :center,
    tcolor
    #theme_palette(:default)[G.groups[i]]
    ; offset=0, textfunc=_short_name,
    #fontargs=(;rotation=-rand(-25:25))
    )
end
# show other labels
plot!()
##
savefig("figures/modularity-tofrom-cc.pdf")
## Try repliedto network
G = _read_final_with_products("fauci-email-repliedto.json")
G = (G..., groups=G.products.simple.modularity, xy=G.products.simple.xy)
drawgraph(G, size=(350,350))
bc = igraph_betweenness(G.A)
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=0.5*sqrt.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1, hover=G.names[dp])

## Test this with a drawing
drawgraph(G, linecolor=:black, size=(550,550))
# show groups ...
using LazySets
for gid = 1:maximum(G.groups)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.75) for i in findall(g->g==gid,G.groups)])
  plot!(h,alpha=0.1, color=gid, linewidth=0)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.5) for i in findall(g->g==gid,G.groups)])
  plot!(h,alpha=0.1, color=gid, linewidth=0)
end
# permute by bc size...
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1)
for i=dp[1:6]
  #annotate!(G.xy[i,1],G.xy[i,2], G.names[i])
  tcolor = G.groups[i] > 16 ? :black : _darker(theme_palette(:default)[G.groups[i]], 0.5)
  showlabel!(G,G.names[i], 9, :center,
    tcolor
    #theme_palette(:default)[G.groups[i]]
    ; offset=0, textfunc=_short_name,
    #fontargs=(;rotation=-rand(-25:25))
    )
end
# show other labels
plot!()
##
savefig("figures/modularity-repliedto.pdf")

##
G = _read_final_with_products("fauci-email-hypergraph-projection.json")
G = (G..., groups=G.products.weighted.modularity, xy=G.products.weighted.xy)

## Test this with a drawing
G = _read_final_with_products("fauci-email-hypergraph-projection.json")
G = _simple_graph(G)
G = (G..., groups=G.products.simple.modularity, xy=G.products.simple.xy)
# drop Fauci
@assert(G.names[1] == "fauci, anthony")
G = (A=G.A[2:end,2:end],groups=G.groups[2:end],names=G.names[2:end],orgs=G.orgs[2:end], xy=G.xy[2:end,:])
bc = igraph_betweenness(G.A)
## try and find better coords
function _fisheye_expand(G)
  # doing this was suggested by Tamara Muzner I believe.
  center = sum(G.xy,dims=1)/size(G.xy,1)
  xy = G.xy .- center
  # convert to polar
  r = sum(xy.^2, dims=2)
  t = atan.(xy[:,2],xy[:,1])
  r = log.(r.+1)
  xy = [r.*cos.(t) r.*sin.(t)]
  (G..., xy=xy)
end

##
G = _fisheye_expand(G)
##
drawgraph(G, linecolor=:black, size=(550,550))
##
# show groups ...
using LazySets
for gid = 1:maximum(G.groups)
  h = ConvexHullArray([Ball2(vec(G.xy[i,:]), 0.25) for i in findall(g->g==gid,G.groups)])
  if gid > 16
    plot!(h,alpha=0.1, color=:black, linewidth=0)
  else
    plot!(h,alpha=0.1, color=gid, linewidth=0)
  end
end
# permute by bc size...
dp = sortperm(bc, rev=true)
scatter!(G.xy[dp,1],G.xy[dp,2],
  markersize=log.(bc[dp].+1).+2, color=G.groups[dp],
  markerstrokecolor=:white,
  markerstrokewidth=1)
for i=dp[1:16]
  #annotate!(G.xy[i,1],G.xy[i,2], G.names[i])
  tcolor = G.groups[i] > 16 ? :black : _darker(theme_palette(:default)[G.groups[i]], 0.5)
  showlabel!(G,G.names[i], 9, :center,
    tcolor
    #theme_palette(:default)[G.groups[i]]
    ; offset=0, textfunc=_short_name,
    #fontargs=(;rotation=-rand(-25:25))
    )
end
# show other labels
plot!()
##
savefig("figures/modularity-hyperproj.pdf")

## Make a table of networks, and the top BC nodes and partitions
graphs = [#"fauci-email-repliedto.json" => "\\texttt{repliedto}",
  #"fauci-email-hypergraph-projection.json" => "\\texttt{hypergraph-projection} without CC",
  "fauci-email-tofrom-5.json"  => "\\texttt{tofrom} without CC",
  "fauci-email-tofrom-cc-5.json"  => "\\texttt{tofrom} with CC"]
topk = 25
results = Dict(map( g-> begin
  G = _read_final_with_products(g)
    G = (G..., groups=G.products.simple.modularity)
  # remove Fauci for comparison...
  if "fauci, anthony" in G.names
    fid = nodeid(G,"fauci, anthony")
    subset = setdiff(1:length(G.names), [fid])
    G = (A=G.A[subset,subset], names=G.names[subset], groups=G.groups[subset])
  end
  G = _simple_graph(G)
  x = igraph_betweenness(G.A)
  p = sortperm(x, rev=true)
  g=>(topk = (G.names[p[1:topk]] .=> G.groups[p[1:topk]]),
    x=x, groups=G.groups, names=G.names)
end, first.(graphs)))

function _namecolor(meta,x)
  id = nodeid(meta, x)
  g = meta.groups[id]
  c = g > 16 ? colorant"black" : theme_palette(:default)[g]
  c = _darker(c, 0.75)
  return "\\textcolor[rgb]{$(round(c.r,digits=2)),$(round(c.g,digits=2)),$(round(c.b,digits=2))}{$(meta.names[id]), $g}"
end

_write_score_table(results,graphs;writescore=false,
  nameformat=_namecolor)
