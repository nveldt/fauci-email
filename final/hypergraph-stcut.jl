##
include("../methods.jl")
include("../methods_hypergraph.jl")
using Measures
H = _read_final_hypergraph("fauci-email-hypergraph.json")
## Project and layout the big hypergraph
xyG = igraph_layout(project_hypergraph(H);random=false)
## Remove big hyperedges
function _make_hypergraph_stcut(H, s::AbstractString, t::AbstractString;maxdeg=10,xy=nothing,
    distribute_hyperedge::Bool=false)
  filt = vec(sum(H.H,dims=2) .<= maxdeg)
  Hf = (H..., H=H.H[filt,:], weights=H.weights[filt], elist=H.elist[findall(filt)])
  ##
  A,ccfilt = project_hypergraph(Hf;distribute_hyperedge) |> A->largest_component(A)
  if xy !== nothing
    xy = xy[ccfilt,:]
  else
    xy = igraph_layout(A; random=false)
  end
  G = (A = A, names = Hf.names[ccfilt], xy = xy,orgs = Hf.orgs[ccfilt])



  Shyp = hyperstcut(Hf.H, nodeid(Hf,s),nodeid(Hf,t),Hf.weights)
  # need to map Shyp to graph via ccfilt ...
  # Shyp =
  @show size(H.H)
  Shyp = begin x=zeros(Int,size(H.H,2)); x[ccfilt.==1] .= 1:sum(ccfilt); x[Shyp] end
  Sgra = stcut(G.A, nodeid(G,s),nodeid(G,t); smallside=false)

  matdraw = triu(G.A,1)
  #lw = 0.5*nonzeros(matdraw)'/(sum(matdraw)/length(nonzeros(matdraw)))
  #lw = 0.5*sqrt.(nonzeros(matdraw)')
  lw = min.(0.5*(nonzeros(matdraw)'), 5)
  plot(draw_graph_lines_vector(matdraw, G.xy;shorten=0.0)..., size=(350,350),
    framestyle=:none, margin=-20mm, label="", alpha=0.5, color=:grey,
    linewidth=lw)

  #=
  Ac = copy(G.A)
  Ac[S,S] .= 0 # remove all edges within S
  Ac[setdiff(1:size(G.A,1),S),setdiff(1:size(G.A,1),S)] .= 0 # edges within Sbar
  dropzeros!(Ac)
  _print_cut(G, S)
  # crude glow simulation
  plot!(draw_graph_lines_vector(Ac, G.xy; shorten=0.01)...,
    linewidth=0.5, color=1, alpha=0.75, label="")
  plot!(draw_graph_lines_vector(Ac, G.xy; shorten=0.02)...,
    linewidth=1, color=1, alpha=0.5, label="")
  plot!(draw_graph_lines_vector(Ac, G.xy; shorten=0.05)...,
      linewidth=2, color=1, alpha=0.25, label="")
  plot!(draw_graph_lines_vector(Ac, G.xy; shorten=0.1)...,
      linewidth=2.5, color=1, alpha=0.1, label="")
  drawset!(G, S,
    marker=:square,
    markersize=2, color=4, markerstrokecolor=:black,  markerstrokewidth=0.5, label="Collins-Conrad")
  drawset!(G, setdiff(1:size(G.A,1),S),
    marker=:circle,
    markersize=2.5, color=2, markerstrokecolor=:black, markerstrokewidth=0.5, label="")
  =#
  JustHype = setdiff(Shyp,Sgra)
  JustGraph = setdiff(Sgra,Shyp)
  both = intersect(Sgra,Shyp)
  rest = setdiff(1:size(G.A,1),[JustHype;JustGraph;both])

  drawset!(G, both,
    marker=:dot,
    markersize=5, color=1, markerstrokecolor=:white,
    markerstrokewidth=0.5, label="")

  drawset!(G, JustHype,
    marker=:square,
    markersize=5, color=2, markerstrokecolor=:white,
    markerstrokewidth=0.5, label="")

  drawset!(G, rest,
    marker=:circle,
    markersize=4, color=:black, markerstrokecolor=:white,
    markerstrokewidth=0.5, label="")

  drawset!(G, JustGraph,
    marker=:downtri,
    markersize=5, color=3, markerstrokecolor=:white,
    markerstrokewidth=0.5, label="")

  showlabel!(G, s, 9, :right; fontargs=(;rotation=-20), offset=2)
  showlabel!(G, t, 9, :left;  fontargs=(;rotation=-10), offset=1)
  showlabel!(G, "redfield",9, :left;  fontargs=(;rotation=-20), offset=1)
  plot!(legend=:bottomleft, legendfontsize=9)
end
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=10, xy=xyG)
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=35, xy=xyG)

##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=5, xy=xyG)
savefig("figures/hypergraph-stcut-maxdeg-5.pdf")
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=10, xy=xyG)
savefig("figures/hypergraph-stcut-maxdeg-10.pdf")
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=15, xy=xyG)
savefig("figures/hypergraph-stcut-maxdeg-15.pdf")
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=25, xy=xyG)
savefig("figures/hypergraph-stcut-maxdeg-25.pdf")
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=100, xy=xyG, distribute_hyperedge=true)
##
_make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=100, xy=xyG, distribute_hyperedge=false)

##
for deg=[5,10,15,25,100]
  _make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=deg, xy=xyG, distribute_hyperedge=false)
  savefig("figures/hypergraph-stcut-maxdeg-$deg.pdf")
end
##
for deg=[5,10,15,25,100]
  _make_hypergraph_stcut(H, "collins", "conrad"; maxdeg=deg, xy=xyG, distribute_hyperedge=true)
  savefig("figures/hypergraph-stcut-distributed-maxdeg-$deg.pdf")
end
