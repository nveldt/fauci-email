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
function _rank_in_others(name,results,keyorder)
  # name - the person to get the rank of
  # results - the dictionary of results
  # myresult - the tag for my result in the results dictionary
  map(key -> begin
      r = results[key]
      if !(name in r.names)
        return (key => missing)
      end
      p = sortperm(r.x, rev=true)
      return (key => findfirst(r.names[p] .== name))
    end, keyorder)
end
function _write_score_table(results, order_and_titles)
  nresults = length(order_and_titles)
  for (key,title) in order_and_titles
    r = results[key]
    println("%")
    println("% -- ", key)
    println("%")
    println("\\begin{tabular}{*{$nresults}{p{16pt}@{}}p{72pt}}")
    println("\\toprule")
    println("\\multicolumn{$(nresults+1)}{c}{$title} \\\\")
    println("\\midrule")
    for (n,v) in r.topk # name and value
      in_others = _rank_in_others(n,results,first.(order_and_titles))
      #@show n, v
      #@show in_others
      for (key_other, rank_in_other) in in_others
        if key_other == key
          print("\\textcolor{LightGray}{")
        end
        print(ismissing(rank_in_other) ? "--" :  rank_in_other)
        if key_other == key
          print("}")
        end
        print(" & ")
      end
      println(n, " \\\\")
    end
    println("\\bottomrule")
    println("\\end{tabular}")
  end
end
_write_score_table(results,graphs)
