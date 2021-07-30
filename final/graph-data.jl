##
include("../methods.jl")
##
using Statistics
#=
the graph data we get is in terms of the simple graph and the weighted graph

simple
- nedges
- max degree
- median degree
- mean degree
- lam2

weighted
- nedges
- nloops
- total volume
- loop volume
- max weighted degree
- median weighted degree
- mean weighted degree
- lam2
=#
function _degree_info(d)
  return (maxdegree=maximum(d), meandegree=mean(d), mediandegree=median(d))
end
function _weighted_data(A)
  @assert all(x->!iszero(x), A.nzval)
  return (
    nverts = size(A,1),
    nedges = Int((nnz(A)-nnz(diag(A)))/2),
    nloops = nnz(diag(A)),
    vol = sum(A),
    loopvol = sum(diag(A)),
    _degree_info(vec(sum(A;dims=1)))...,  # read this out into the return val
    lam2 = fiedler_vector(A)[2]
  )
end
function _simple_data(A)
  @assert nnz(diag(A)) == 0
  @assert all(x->x==1, A.nzval)
  return _weighted_data(A) # includes redundant data, ohwell...
end
##
function _graph_data(G)
  A = G.A
  B = spones!(dropzeros!(A - Diagonal(A)))
  # make simple graph
  nverts = size(A,1)
  simple = _simple_data(B)
  weighted = _weighted_data(A)
  return (nverts=nverts, simple=simple, weighted=weighted)
end
function _write_graph_data_table_latex(G, gname)
  data = _graph_data(G)
  s = data.simple
  w = data.weighted

  str = join([#"",
    string(gname),
    string(data.nverts),
    string(s.nedges),
    string(Int(s.maxdegree)),
    string(round(s.meandegree,digits=1)),
    string(Int(s.mediandegree)),
    string(round(s.lam2,digits=4)),
    #"",
    #string(w.nedges),
    string(w.nloops),
    string(Int(w.vol)),
    string(Int(w.loopvol)),
    string(Int(w.maxdegree)),
    string(round(w.meandegree,digits=1)),
    string(Int(w.mediandegree)),
    string(round(w.lam2,digits=4))*"\\\\"],
    " & ")
  println(str)
end
function _write_graph_data_table_markdown(G, gname)
  data = _graph_data(G)
  s = data.simple
  w = data.weighted

  str = join(["",
    string(gname),
    string(data.nverts),
    string(s.nedges),
    string(Int(s.maxdegree)),
    string(round(s.meandegree,digits=1)),
    string(Int(s.mediandegree)),
    string(round(s.lam2,digits=4)),
    "",
    #string(w.nedges),
    string(w.nloops),
    string(Int(w.vol)),
    string(Int(w.loopvol)),
    string(Int(w.maxdegree)),
    string(round(w.meandegree,digits=1)),
    string(Int(w.mediandegree)),
    string(round(w.lam2,digits=4)),
    ""],
    " | ")
  println(str)
end
graphs = ["fauci-email-repliedto.json" => "\\texttt{repliedto-nofauci}",
  "fauci-email-hypergraph-projection.json" => "\\texttt{hypergraph-proj} w/o CC",
  "fauci-email-hypergraph-projection-cc.json" => "\\texttt{hypergraph-proj} with CC",
  "fauci-email-tofrom-5.json"  => "\\texttt{tofrom-nofauci} w/o CC",
  "fauci-email-tofrom-cc-5.json"  => "\\texttt{tofrom-nofauci} with CC",]
println()
map(x -> begin
  G = _read_final(x[1])
  _write_graph_data_table_latex(G, x[2])
end, graphs)

##results = Dict(map( g-> begin
println()
map(x -> begin
  G = _read_final(x[1])
  _write_graph_data_table_markdown(G, x[2])
end, graphs)
