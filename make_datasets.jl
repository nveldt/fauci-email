include("methods.jl")

function _write_list(io, listname::String, items; last::Bool=false, map::Function=string)
  nitems = length(items)
  write(io, "\"$listname\": [", "\n")
  for (ii,i) in enumerate(items)
    write(io, map(i), ii < nitems ? ",\n" : "\n")
  end
  write(io, last ? "]\n" : "],\n") # skip comma if last...
end
function _write_edgedata(io, A::SparseMatrixCSC, nedges::Integer; last::Bool=false)
  write(io, "\"edgedata\": [", "\n")
  # ugh, annoying to get this precise for JSON with the last comma missing...
  edata = zip(findnz(A)...)
  @assert(length(edata) == nedges)
  for (ei,(i,j,w)) in enumerate(edata)
    write(io, string(i-1), ", ", string(j-1), ", ", string(w), ei < nedges ? ",\n" : "\n")
  end
  write(io, last ? "]\n" : "],\n") # skip comma if last...
end


""" Write a custom json output method to make it easy to parse these
files without depending on a full json library... but also make them
JSON so if you want to use that, it's easy! The format is:
```
{
  "vertices": <number of vertices>,
  "edges": <number of edges>,
  "edgedata": [
    <src1>, <dst1>, <weight1>,
    <src2>, <dst2>, <weight2>,
    ...
    <src_number_of_edges>, <dst_number_of_edges>, <weight_number_of_edges>
  ],
  "labels": [
    <list of labels, one per vertex>
  ],
  "orgs": [
    <list of organizations, one per vertex"
  ]
}
```
so that you could parse these as follows (in Python, say...)
with open("gdata.json", "r") as f:
  f.readline() # read the first '{'
  nverts = int(f.readline().split(':')[1].split(',')[0])
  nedges = int(f.readline().split(':')[1].split(',')[0])
  f.readline() # read "edgedata"
  src, dst, weights = [],[],[]
  for _ in range(nedges):
    einfo = f.readline().split(",")
    src.append(int(einfo[0]))
    dst.append(int(einfo[1]))
    weights.append(int(einfo[2]))
  f.readline() # read end array
  f.readline() # read label array start
  labels = []
  for _ in range(nverts)
    labels.append(f.readline().strip().strip('"'))
  orgs = []
  for _ in range(nverts)
    orgs.append(int(f.readline()))
```
all information is zero-indexed...
an undirected edge is repeated twice, one for each direction.
If you need to translate into a SNAP graph and drop weights, say, this is easy with
shell tools
```
\$ tail -n +5 fauci-email-tofrom-5.json | sed -n '/],/q;p' | sed 's/,//g' | cut -f1,2 -d" " | less
```
"""
function _write_simple_json(filename::String, G::NamedTuple)
  nverts = size(G.A,1)
  nedges = nnz(G.A)
  labels = G.names
  orgs = G.orgs
  open(filename, "w") do f
    write(f, "{\n")
    write(f, "\"vertices\": ", string(nverts), ",\n")
    write(f, "\"edges\": ", string(nedges), ",\n")
    _write_edgedata(f, G.A, nedges)
    @assert(length(labels) == nverts)
    _write_list(f, "labels", labels; map=JSON.json)
    @assert(length(orgs) == nverts)
    _write_list(f, "orgs", orgs; last=true)
    write(f, "}\n")
  end
end
##
G = _build_email_tofrom_graph(data; maxset=5, keepfauci=false)
_write_simple_json("fauci-email-graph-tofrom-nofauci-nocc-5.json", G)

##
G = _build_email_tofrom_graph(data; maxset=5, keepfauci=false, keepcc=true)
_write_simple_json("fauci-email-graph-tofrom-nofauci-cc-5.json", G)
##
G = _build_email_repliedto_graph(data; keepfauci=false)
_write_simple_json("fauci-email-graph-repliedto-nofauci.json", G)
##
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients"), mindegree=2)
_write_simple_json("fauci-email-graph-hypergraph-projection-nocc.json", G)
##
G = _build_email_hypergraph_projection(data;
  hyperedgeparts=("sender","recipients","cc"), mindegree=2)
_write_simple_json("fauci-email-graph-hypergraph-projection-cc.json", G)
##
Tcc_ids = temporal_reachability(data) |> R->min.(R.R, R.R') |> simple_clique_heuristic

T = build_temporal_graphs(data; subset=Tcc_ids)
function _write_json_graph_sequence(filename::String, G::NamedTuple)
  nverts = size(G.T[1][2],1)
  ngraphs = length(G.T)
  labels = G.names
  orgs = G.orgs
  open(filename, "w") do f
    write(f, "{\n")
    write(f, "\"vertices\": ", string(nverts), ",\n")
    write(f, "\"graphs\": ", string(length(G.T)), ",\n")
    write(f, "\"sequencedata\": [")
    for (gi, (date, g)) in enumerate(G.T)
      nedges = nnz(g)
      write(f, "{\n")
      write(f, "\"edges\": ", string(nedges), ",\n")
      _write_edgedata(f, g, nedges; last=true)
      write(f, gi < length(G.T) ? "},\n" : "}\n")
    end

    write(f, "],", "\n")
    @assert issorted(first.(G.T))
    _write_list(f, "dates", first.(G.T); map=JSON.json)
    @assert(length(labels) == nverts)
    _write_list(f, "labels", labels; map=JSON.json)
    @assert(length(orgs) == nverts)
    _write_list(f, "orgs", orgs; last=true)
    write(f, "}\n")
  end
end
_write_json_graph_sequence("fauci-email-temporalgraph-tofrom.json", T)

## Construct a hypergraph and save it.
include("methods_hypergraph.jl")
kf = false
ms = 100
mindeg = 5
parts=("sender","recipients")
H = _build_email_hypergraph(data;hyperedgeparts=parts,maxset=ms, keepfauci=kf,mindegree = mindeg)

##
function _write_hyperedgedata(io, A::SparseMatrixCSC, nedges::Integer; last::Bool=false)
  write(io, "\"hyperedgedata\": [", "\n")
  # ugh, annoying to get this precise for JSON with the last comma missing...
  edata = zip(findnz(A)...)
  @assert(length(edata) == nedges)
  for (ei,(i,j,w)) in enumerate(edata)
    write(io, string(j-1), ", ", string(i-1), ei < nedges ? ",\n" : "\n")
  end
  write(io, last ? "]\n" : "],\n") # skip comma if last...
end
function _write_hypergraph_json(filename::String, G::NamedTuple)
  nverts = size(G.H,2)
  nedges = size(G.H,1)
  nincidence = nnz(G.H)

  labels = G.names
  orgs = G.orgs
  weights = G.weights
  open(filename, "w") do f
    write(f, "{\n")
    write(f, "\"vertices\": ", string(nverts), ",\n")
    write(f, "\"hyperedges\": ", string(nedges), ",\n")
    write(f, "\"incidences\": ", string(nincidence), ",\n")
    _write_hyperedgedata(f, copy(G.H'), nincidence)
    @assert(length(weights) == nedges)
    _write_list(f, "weights", weights)
    @assert(length(labels) == nverts)
    _write_list(f, "labels", labels; map=JSON.json)
    @assert(length(orgs) == nverts)
    _write_list(f, "orgs", orgs; last=true)
    write(f, "}\n")
  end
end
_write_hypergraph_json("fauci-email-hypergraph.json", H)
