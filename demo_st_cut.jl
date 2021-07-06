include("methods.jl")

## Nate's large set analysis
G = _build_email_tofrom_graph(data; maxset=5, keepfauci=false) |> igraph_layout
drawgraph(G)
drawset!(G, stcut(G, "collins", "conrad"))


## David's media set analysis
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"), mindegree=2) |> igraph_layout
drawgraph(G)
drawset!(G, stcut(G, "fauci", "conrad"))