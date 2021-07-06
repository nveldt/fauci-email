include("methods.jl")

## Nate's large set analysis
G = _build_email_tofrom_graph(data; maxset=5, keepfauci=false) |> igraph_layout
drawgraph(G)
# This syntax flourish saves S1N, the names in the set, and gives S1 to drawset!
drawset!(G, begin S1 = stcut(G, "collins", "conrad");
                  S1N = G.names[setdiff(1:length(G.names), S1)]; # need to flip set to analysze with others below...
                  S1; end )

## ... without the maxset restriction ... but with weights
G = _build_email_tofrom_graph(data; keepfauci=false, emailweight=true) |> igraph_layout
drawgraph(G)
drawset!(G, begin S2 = stcut(G, "collins", "conrad"); S2N = G.names[S2]; S2; end)

## ... in the hypergraph projection  ...
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"),
  keepfauci = false, emailweight = true) |> igraph_layout
drawgraph(G)
drawset!(G, begin S3 = stcut(G, "collins", "conrad"); S3N = G.names[S3]; S3; end)

## ... compare the sets ...
function compare_sets(args...)
  setdiffs = zeros(length(args),length(args))
  for i=1:length(args)
    for j=1:length(args)
      if i != j
        setdiffs[i,j] = length(setdiff(args[i], args[j]))
      end
    end
  end
  return setdiffs
end
compare_sets(S1N,S2N,S3N)

## David's media set analysis
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"), mindegree=2) |> igraph_layout
drawgraph(G)
drawset!(G, stcut(G, "fauci", "conrad"))

## ... in the tofrom graph ...
G = _build_email_tofrom_graph(data, mindegree=2) |> igraph_layout
drawgraph(G)
drawset!(G, stcut(G, "fauci", "conrad"))

## And this is also present in the full graphs
G1 = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients")) |> igraph_layout
FC1 = stcut(G1, "fauci", "conrad") |> x-> G1.names[x]

G2 = _build_email_tofrom_graph(data) |> igraph_layout
drawset!(G2, stcut(G2, "fauci", "conrad"))
FC2 = stcut(G2, "fauci", "conrad") |> x-> G2.names[x]

#  ... and this is the difference
@show setdiff(FC1,FC2)
@show setdiff(FC2,FC1)
