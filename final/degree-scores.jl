# to run,
#
include("../methods.jl")
##
graphs = ["fauci-email-repliedto.json" => "\\texttt{repliedto}",
  "fauci-email-hypergraph-projection.json" => "\\texttt{hypergraph-projection} without CC",
  "fauci-email-hypergraph-projection-cc.json" => "\\texttt{hypergraph-projection} with CC",
  "fauci-email-tofrom-5.json"  => "\\texttt{tofrom} without CC",
  "fauci-email-tofrom-cc-5.json"  => "\\texttt{tofrom} with CC",]
topk = 10
results = Dict(map( g-> begin
  G = _read_final(g)
  # remove Fauci for comparison...
  if "fauci, anthony" in G.names
    fid = nodeid(G,"fauci, anthony")
    subset = setdiff(1:length(G.names), [fid])
    G = (A=G.A[subset,subset], names=G.names[subset])
  end
  x = vec(sum(G.A;dims=1))
  p = sortperm(x, rev=true)
  g=>(topk = (G.names[p[1:topk]] .=> x[p[1:topk]]), x=x, names=G.names)
end, first.(graphs)))

##
map( g-> begin
  G = _read_final(g)
  @show g, issymmetric(G.A)
end, first.(graphs))

##
_write_score_table(results,graphs)
