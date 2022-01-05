# to run,
#
include("methods.jl")
##
graphs = ["fauci-email-graph-repliedto-nofauci.json" => "\\texttt{repliedto-nofauci}",
  "fauci-email-graph-hypergraph-projection-nocc.json" => "\\texttt{hypergraph-proj} w/o CC",
  "fauci-email-graph-hypergraph-projection-cc.json" => "\\texttt{hypergraph-proj} with CC",
  "fauci-email-graph-tofrom-nofauci-nocc-5.json"  => "\\texttt{tofrom-nofauci} w/o CC",
  "fauci-email-graph-tofrom-nofauci-cc-5.json"  => "\\texttt{tofrom-nofauci} with CC"]
topk = 10
results = Dict(map( g-> begin
  G = _read_final(g)
  # remove Fauci for comparison...
  if "fauci, anthony" in G.names
    fid = nodeid(G,"fauci, anthony")
    subset = setdiff(1:length(G.names), [fid])
    G = (A=G.A[subset,subset], names=G.names[subset])
  end
  x = pagerank(G.A, 0.85)
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
