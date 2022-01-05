##
include("../methods.jl")
##
G = _read_final_with_products("fauci-email-graph-tofrom-nofauci-nocc-5.json")
using Measures
##
function _print_cut(G, S)
  Ac = copy(G.A)
  Ac[S,S] .= 0 # remove all edges within S
  Ac[setdiff(1:size(G.A,1),S),setdiff(1:size(G.A,1),S)] .= 0 # edges within Sbar
  dropzeros!(Ac)
  cutedges = zip(findnz(triu(Ac,1))[1:2]...)
  for e in cutedges
    println("$(G.names[e[1]]) & $(G.names[e[2]]) \\\\")
  end
end
function _st_plot(G, s::AbstractString, t::AbstractString)
  #drawgraph(G, pointcolor=:none, label="", size=(350,350))
  @show extrema(G.A.nzval)
  plot(draw_graph_lines(G.A, G.xy)..., size=(350,350),
    framestyle=:none, margin=-20mm, label="", alpha=0.5, color=:black,
    linewidth=0.75)
  S = stcut(G, s, t)
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
  showlabel!(G, s, 7, :left; rotation=15, offset=2)
  showlabel!(G, t, 7, :left; rotation=15, offset=1)
  plot!(legend=:bottomleft, legendfontsize=9)
end
_st_plot(G, "collins", "conrad")
showlabel!(G,"kadlec",7)
showlabel!(G,"lane",7)
showlabel!(G,"marston",7)
showlabel!(G,"stover",7)

##
_st_plot(G, "collins", "conrad")
plot!(legend=:none)
##
savefig("figures/st-cut-tofrom-collins-conrad-modularity.pdf")

##
_st_plot(igraph_layout(G;random=false), "collins", "conrad")
plot!(legend=:none)
##
savefig("figures/st-cut-tofrom-collins-conrad-force.pdf")

## make the graph simple
function _st_plot_hyper(G, s::AbstractString, t::AbstractString)
  #drawgraph(G, pointcolor=:none, label="", size=(350,350))
  plot(draw_graph_lines(G.A, G.xy)..., size=(350,350),
    framestyle=:none, margin=-20mm, label="", alpha=0.1, color=:black, linewidth=0.5)
  S = stcut(G, s, t)
  Ac = copy(G.A)
  Ac[S,S] .= 0 # remove all edges within S
  Ac[setdiff(1:size(G.A,1),S),setdiff(1:size(G.A,1),S)] .= 0 # edges within Sbar
  triu!(Ac)
  dropzeros!(Ac)
  AS = copy(G.A)
  AS[setdiff(1:size(G.A,1),S),setdiff(1:size(G.A,1),S)] .= 0 # edges within Sbar
  triu!(AS)
  dropzeros!(AS)
  println.(G.names[S])
  # Conrad edges...
  plot!(draw_graph_lines(AS, G.xy)..., size=(350,350),
    framestyle=:none, label="", alpha=0.5, color=:black, linewidth=0.75)
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
    markersize=2, color=4, markerstrokecolor=:black,  markerstrokewidth=0.5, label="")
  drawset!(G, setdiff(1:size(G.A,1),S),
    marker=:circle,
    markersize=2.5, color=2, markerstrokecolor=:black, markerstrokewidth=0.5, label="")
  showlabel!(G, s, 7, :left, :black; fontargs=(;rotation=35), offset=2)
  showlabel!(G, t, 7, :left, :black; fontargs=(;rotation=35), offset=2)

  plot!(legend=:bottomleft, legendfontsize=9)
end
G = _read_final_with_products("fauci-email-graph-hypergraph-projection-nocc.json")
_st_plot_hyper(G, "conrad", "fauci")
plot!(legend=:none)
##
savefig("figures/st-cut-hyperproj-conrad-fauci-modularity.pdf")
##

##
_st_plot_hyper(igraph_layout(G;random=false), "conrad", "fauci")
plot!(legend=:none)
##
savefig("figures/st-cut-hyperproj-conrad-fauci-force.pdf")
