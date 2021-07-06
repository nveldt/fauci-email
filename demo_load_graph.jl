include("methods.jl")

## To-from graph loading
keepfauci = false
keepcc = false
eweight = false     # scales by 1/(|e|-1), |e| = number people in email
maxset = 5
mindeg = 0
G = _build_email_tofrom_graph(data; maxset=maxset, keepfauci=keepfauci,
    mindegree = mindeg,keepcc = keepcc,emailweight = eweight) |> igraph_layout

drawgraph(G,markersize = 6)    

## Hypergraph projection graph
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"), mindegree=2,
  keepfauci = keepfauci, emailweight = eweight) |> igraph_layout
drawgraph(G)