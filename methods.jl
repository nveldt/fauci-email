## Data
using LinearAlgebra
using MatrixNetworks
using JSON
using DelimitedFiles
using Printf

data = JSON.parsefile("fauci-email-graph.json")
include("include/FlowSeed.jl")
include("include/PushRelabelMaxFlow.jl")

##
include("methods_temporal.jl")

function stcut(A::SparseMatrixCSC, s::Int, t::Int; smallside::Bool = true)
    F = maxflow(Float64.(A),s,t,0.0)
    S = source_nodes_min(F)
    if smallside
      if length(S) > size(A,1)/2
        S = setdiff(1:size(A,1), S)
      end
    end
    return S
end

function bigsplits(A,minsize::Int=15)
    n = size(A,1)
    bigsets = NTuple{2,Int}[]
    sizes = Int[]
    for i=1:n
        for j=i+1:n
            S = stcut(A,i,j;smallside=true)
            if length(S) > minsize
                push!(bigsets, tuple(i,j))
                push!(sizes, length(S))
            end

        end
    end
    # reshape to matrix
    seeds = zeros(length(bigsets), 2)
    for i=1:length(bigsets)
        seeds[i,1] = bigsets[i][1]
        seeds[i,2] = bigsets[i][2]
    end
    return seeds, sizes
end


##
using PyCall
using Conda
using Random
const igraph = pyimport_conda("igraph","python-igraph","conda-forge")
const pyrandom = pyimport("random")
function igraph_layout(A::SparseMatrixCSC{T}, layoutname::AbstractString="fr";
    random::Bool=true) where T
    ei,ej,ew = findnz(A)
    edgelist = [(ei[i]-1,ej[i]-1) for i = 1:length(ei)]
    nverts = size(A)
    G = igraph.Graph(nverts, edges=edgelist, directed=false)
    if random
      xy = G.layout(layoutname)
    else
      pyrngstate = pyrandom.getstate()
      pyrandom.seed(0)
      xy = G.layout(layoutname)
      pyrandom.setstate(pyrngstate)
    end
    xy = [Float64(xy[i][j]) for i in 1:length(xy),  j in 1:length(xy[1])]
end

function _mindegree_and_cc_filter(A;mindegree=0,weightdegree=false)
  if weightdegree
    d = vec(sum(A;dims=2))
  else
    Aones = copy(A)
    fill!(Aones.nzval, 1)
    d = vec(sum(Aones;dims=2))
  end
  filt = d .>= mindegree
  Ad = A[filt,filt]
  filtids = findall(filt) # get the ids
  Acc,ccfilt = largest_component(Ad)
  return dropzeros!(Acc), filtids[ccfilt]
end

function subset_and_mindegree_cc_filter(G::NamedTuple, subset;
    mindegree=1, connected::Bool=true, weightdegree=false)

  As = G.A[subset,subset]
  names = G.names[subset]
  orgs = G.orgs[subset]

  if weightdegree
    d = vec(sum(As;dims=2))
  else
    Aones = copy(As)
    fill!(Aones.nzval, 1)
    d = vec(sum(Aones;dims=2))
  end

  filt = d .>= mindegree
  Ad = As[filt,filt]
  filtids = findall(filt) # get the ids
  Acc,ccfilt = largest_component(Ad)

  return (G..., A=Acc, names=names[filtids[ccfilt]], orgs=orgs[filtids[ccfilt]])
end

function _build_email_hypergraph_projection(data;
    maxset::Int=typemax(Int), mindegree::Int=0, keepfauci=true,
    hyperedgeparts=("sender","recipients","cc"),
    emailweight::Bool=false)
  edges = Tuple{Int,Int}[]
  weights = Float64[]
  emails=data["emails"]
  names=data["names"]
  idfauci=findfirst(names .== "fauci, anthony")-1 # offset
  for thread in emails
    for email in thread
      people = Set{Int}()
      for group in hyperedgeparts
        map(p->push!(people, p), email[group])
      end
      # project the hyperedge
      if length(people) <= maxset
        for (pi,pj) in Iterators.product(people, people)
          if pi > pj
            if keepfauci || (keepfauci == false && pj!=idfauci && pi!=idfauci)
              push!(edges, (pi+1,pj+1))
              push!(weights, 1/binomial(length(people),2))
            end
          end
        end
      end
    end
  end
  A = sparse(first.(edges), last.(edges), emailweight ? weights : 1, length(names), length(names))
  A = dropzeros!(max.(A,A'))
  println("Raw data has ")
  println("  $(size(A,1)) nodes and $(nnz(A)÷2) edges")
  Af,filt = _mindegree_and_cc_filter(A;mindegree)
  println("Degree (mindegree=$mindegree) and CC filtered data has ")
  println("  $(size(Af,1)) nodes and $(nnz(Af)÷2) edges")

  return (A=Af, names=names[filt], orgs = data["clusters"][filt], filt=filt,
      maxset, mindegree, keepfauci, hyperedgeparts, emailweight)
end

#F = _build_email_hypergraph_projection(data;
  #mindegree=2, hyperedgeparts=("sender", "recipients"))
##

# This is also used in methods-temporal.jl too...
function _build_tofrom_edges(emails;
            edges = Tuple{Int,Int}[],
            weights = Float64[],
            keepcc::Bool=false,
            keepfauci::Bool=true,
            idfauci::Int=-1,
            maxset::Int = typemax(Int))
  @assert(keepfauci==true || (keepfauci==false && idfauci != -1))
  for email in emails
    from=email["sender"]
    if keepcc
      others = [email["recipients"]; email["cc"]]
    else
      others = email["recipients"]
    end
    if length(others) <= maxset
      for p in others
        pi = from
        pj = p # weird system to keep same check as below...
        if keepfauci || (keepfauci == false && pj!=idfauci && pi!=idfauci)
          push!(edges, (pi+1,pj+1))
          push!(weights, 1/length(others))
        end
      end
    end
  end
  return edges, weights
end

function _build_email_tofrom_graph(data;
    maxset::Int=typemax(Int), mindegree::Int=0, keepfauci=true,
    emailweight::Bool=false, keepcc=false)
  names = data["names"]
  edges = Tuple{Int,Int}[]
  weights = Float64[]
  emails=data["emails"]
  idfauci=findfirst(data["names"] .== "fauci, anthony")-1
  for thread in emails
    # note order among the named args doesn't matter, thanks Julia devs!!
    _build_tofrom_edges(thread; edges, weights, keepcc, keepfauci, maxset, idfauci)
  end
  A = sparse(first.(edges), last.(edges), emailweight ? weights : 1,
        length(names), length(names))
  A = dropzeros!(max.(A,A'))
  println("Raw data has ")
  println("  $(size(A,1)) nodes and $(nnz(A)÷2) edges")
  Af,filt = _mindegree_and_cc_filter(A;mindegree)
  println("Degree (mindegree=$mindegree) and CC filtered data has ")
  println("  $(size(Af,1)) nodes and $(nnz(Af)÷2) edges")
  return (A=Af, names=names[filt], orgs = data["clusters"][filt], filt=filt,
      maxset, mindegree, keepfauci, emailweight)
end
##
""" Build a replied-to graph. In this graph, two nodes are connected
if x replied to an email from y. This is built by looking at each
  thread in reverse order (so oldest message first). Take the
  sender of message 1... and then the sender of message 2... we put in an edge
  from sender 1 to sender 2 if sender 1 is either a recipient or CC of message 2.
  This analysis is done for all temporally ordered pairs of emails
    in the thread.

  Note that ideally we would run this on the threads without duplicate emails.
"""
function _build_email_repliedto_graph(data;
    mindegree::Int=0, keepfauci=true)
  names = data["names"]
  edges = Tuple{Int,Int}[]
  weights = Float64[]
  emailweight = false
  emails=data["emails"]
  idfauci=findfirst(data["names"] .== "fauci, anthony")-1
  for thread in emails
    nemails = length(thread)
    for eid1=1:nemails # eid1 is the older email
      for eid2=1:eid1-1  #eid2 are all possible emails that reply to eid1
        e1 = thread[eid1]
        e2 = thread[eid2]
        s1 = e1["sender"]
        s2 = e2["sender"]
        # s2 is the newer email... so is s1 a recipeient or CC of s2?
        if s1 in e2["recipients"] || s1 in e2["cc"]
          if keepfauci || (keepfauci == false && s1!=idfauci && s2!=idfauci)
            push!(edges, (s1+1,s2+1))
          end
        end
      end
    end
  end
  A = sparse(first.(edges), last.(edges), emailweight ? weights : 1,
        length(names), length(names))
  A = dropzeros!(max.(A,A'))
  println("Raw data has ")
  println("  $(size(A,1)) nodes and $(nnz(A)÷2) edges")
  Af,filt = _mindegree_and_cc_filter(A;mindegree)
  println("Degree (mindegree=$mindegree) and CC filtered data has ")
  println("  $(size(Af,1)) nodes and $(nnz(Af)÷2) edges")
  return (A=Af, names=names[filt], orgs = data["clusters"][filt], filt=filt,
      mindegree, keepfauci, emailweight)
end
#F = _build_email_tofrom_graph(data; mindegree=0, keepfauci=false, maxset=5)

##
function _simple_graph(G::NamedTuple)
  return (G..., A = spones!(dropzeros!(G.A - Diagonal(G.A))))
end

##
function densify_graph_groups(A,groups; within_group_deg::Int=ceil(Int,length(unique(groups))),rng=Random.MersenneTwister(0))
  A1 = spones(A)
  ngroups = maximum(groups)
  @assert minimum(groups) >= 1
  @assert length(unique(groups)) <= ngroups
  G = sparse(1:length(groups),groups,1,size(A,1),maximum(groups))
  println("density_graph_groups -- within_group_deg=$(within_group_deg)")
  for gi = 1:ngroups
    Gv = findnz(G[:,gi])[1]
    # generate a little local ER graph with within_group_deg
    # current group degree
    gd = sum(A1[Gv,Gv])/length(Gv)
    target_deg = within_group_deg

    if gd >= target_deg
      # move it halfway to max possible
      target_deg = 0.5*gd + 0.5*length(Gv)
    end
    erp = min((target_deg - gd)/length(Gv),1)
    println("group $(gi) size=$(length(Gv))  within_avg_deg = $(gd)   target_within=$(target_deg)  erprob=$(erp)")
    Gr = sprand(rng, length(Gv),length(Gv),erp/2) # we will symmetrize to double erp
    Gr = max.(Gr,Gr')
    Gr = spones!(Gr)
    Gr = Gr - Diagonal(Gr)
    A1[Gv,Gv] .+= Gr # add in new edges...
  end
  return spones!(A1)
end

function group_layout(G::NamedTuple,groups;layout="fr", kwargs...)
  xy = igraph_layout(densify_graph_groups(G.A,groups;kwargs...),layout)
  return (G..., groups=groups, xy=xy)
end

## Work with differences, using names to align...
function _realign_set_to_names(G::NamedTuple, S, names)
  idmap_newnames = Dict{String,Int}(names .=> 1:length(names))
  return map(x->idmap_newnames[x], G.names[S])
end
function _realign_to_names(G::NamedTuple, names)
  idmap_newnames = Dict{String,Int}(names .=> 1:length(names))
  idmap = Dict{Int,Int}( 1:length(G.names) .=> map(x->idmap_newnames[x], G.names) )
  I,J,V = findnz(G.A)
  mapids = x->idmap[x]
  H = sparse(map(mapids,I),map(mapids,J),V,length(names),length(names))
  return (A=H, names=names)
end
""" In general, G1 should be bigger than G2... but this isn't a requirement... """
function graph_difference(G1::NamedTuple, G2::NamedTuple)
  add_to_G1 = setdiff(G2.names,G1.names) # what is in G2 but not G1...
  Dnames = append!(G1.names,add_to_G1) # names for the difference...
  H1 = _realign_to_names(G1,Dnames)
  H2 = _realign_to_names(G2,Dnames)
  if length(add_to_G1) == 0
    # then everything in G2 is in G1... so we can reuse all the metadata...
    @assert Dnames == G1.names
    return (G1...,A=H1.A - H2.A) # replace the graph with the difference..
  else
    return (A=(H1.A - H2.A), names=Dnames)
  end
end

function add_orgs(G::NamedTuple, data)
  # index orgs
  name2org = Dict{String,Int}(data["names"] .=> data["clusters"])
  orgs = map(x->name2org, G.names)
  return (G...,orgs=orgs)
end




## Helpers
nodeid(G,p) = G.names |> x-> startswith.(x, p) |> x->findfirst(x)
function spones!(A::SparseMatrixCSC)
  fill!(A.nzval, 1)
  return A
end
spones(A::SparseMatrixCSC) = spones!(copy(A))
stcut(G::NamedTuple, s::String, t::String) = stcut(G.A, nodeid(G,s), nodeid(G,t))
stcut(f::Function, G::NamedTuple, s::String, t::String) = stcut(f(G.A), nodeid(G,s), nodeid(G,t))

function igraph_layout(G::NamedTuple,layout="fr";random::Bool=true)
  xy = igraph_layout(G.A, layout; random)
  return (G..., xy=xy)
end

using Plots
function drawgraph(G::NamedTuple; pointcolor=:auto,
    shownames=false, namecolor=:black, shorten::Float64=0.0,
    kwargs...)
  p = plot(draw_graph_lines_vector(G.A,G.xy;shorten)...;label=false, alpha=0.1,
      framestyle=:none, linecolor=1,
      size=(800,800), kwargs...)
  if pointcolor==:auto
    if haskey(G,:groups)
      pointcolor=:groups
    elseif haskey(G, :orgs)
      pointcolor=:orgs
    else
      pointcolor=:none
    end
  end
  # G[:symbol] = G.symbol, some fancy syntax...
  pointcolor == :none ?  markercolor=:black : markercolor=G[pointcolor]
  scatter!(p,G.xy[:,1],G.xy[:,2],markercolor=markercolor,markersize=2,
      markerstrokewidth=0, legend=false, hover=G.names; kwargs...)
  if shownames
    for i=1:length(G.names)
      annotate!(p, G.xy[i,1],G.xy[i,2], (split(G.names[i], ",")[1], 7, namecolor))
    end
  end
  return p
end

function draw_graph_lines(A::SparseMatrixCSC, xy)
  if issymmetric(A)
      ei,ej = findnz(triu(A,1))[1:2]
  else
      ei,ej = findnz(A)[1:2]
  end
  # find the line segments
  lx = zeros(0)
  ly = zeros(0)
  for nz=1:length(ei)
      src = ei[nz]
      dst = ej[nz]
      push!(lx, xy[src,1])
      push!(lx, xy[dst,1])
      push!(lx, Inf)
      push!(ly, xy[src,2])
      push!(ly, xy[dst,2])
      push!(ly, Inf)
  end
  return lx, ly
end

# We need this one to show different linewidths for each line <shrug>
function draw_graph_lines_vector(A::SparseMatrixCSC, xy; shorten::Float64=0)
  ei,ej = findnz(A)[1:2]
  # find the line segments
  lx = Vector{Float64}[]
  ly = Vector{Float64}[]
  for nz=1:length(ei)
    src = ei[nz]
    dst = ej[nz]
    if shorten == 0
      push!(lx, [xy[src,1], xy[dst,1]])
      push!(ly, [xy[src,2], xy[dst,2]])
    else
      # view xy as a line, then
      a = shorten
      b = 1-shorten
      push!(lx, [b*xy[src,1]+a*xy[dst,1], a*xy[src,1]+b*xy[dst,1]])
      push!(ly, [b*xy[src,2]+a*xy[dst,2], a*xy[src,2]+b*xy[dst,2]])
    end
  end
  return lx, ly
end

function drawset!(G::NamedTuple, S; kwargs...)
  scatter!(G.xy[S,1],G.xy[S,2];hover=G.names[S], markerstrokewidth=0,kwargs...)
end

function showlabel!(G::NamedTuple, name::String, args...; offset::Int=0, fontargs=(;), kwargs...)
  id = nodeid(G,name)
  annotate!(G.xy[id,1],G.xy[id,2],
    Plots.text(repeat(" ",offset)*G.names[id], args...; fontargs...); kwargs...)
end

##
function print_email(email, names; text::Bool=true)
  fromname = names[email["sender"]+1]
  tonames = join(map(x-> names[x+1], email["recipients"]), "; ")
  ccnames = join(map(x-> names[x+1], email["cc"]), "; ")
  println("From: $fromname")
  println("To: $tonames")
  println("CC: $ccnames")
  println("Date: ", email["time"])
  println("Subject: ", email["subject"])
  if text
    println("Text: ", email["text"])
  end
end
function print_email(data::Dict, id::Tuple{Int,Int})
  print_email(data["emails"][id[1]][id[2]], data["names"])
end
function print_thread(thread, names::Vector)
  print_email(thread[1], names)
  for email in thread[2:end]
    println(" *in reply to* ")
    print_email(email, names)
  end
end
function print_thread(data::Dict, tid::Int)
  print_thread(data["emails"][tid], data["names"])
end

## ranking display helpers
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
    println("\\begin{tabular}{*{$nresults}{p{16pt}@{}}p{112pt}}")
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
      println("$n, $(round(v,digits=6))", " \\\\")
    end
    println("\\bottomrule")
    println("\\end{tabular}")
  end
end


##
function _edgedata_to_sparse(gdata, n::Integer)
  if haskey(gdata, "edges") && haskey(gdata, "edgedata")
    m = gdata["edges"]
    X = reshape(Float64.(gdata["edgedata"]), 3, m)
    A = sparse(Int.(X[1,:]).+1, Int.(X[2,:]).+1, X[3,:], n, n)
  else
    throw(ArgumentError("dictionary doesn't have the right keys"))
  end
end
function _read_final(fn::AbstractString)
  gdata = JSON.parsefile(fn)
  n = gdata["vertices"]
  A = _edgedata_to_sparse(gdata, n)
  # find all instances of edgedata and convert to adjaency matrices...
  return (A=A, names=string.(gdata["labels"]), orgs=Int.(gdata["orgs"]))
end

function _read_products(fn::AbstractString)
  p = JSON.parsefile(fn)
  if haskey(p,"xy")
    p["xy"] = hcat(p["xy"]...)
  end
  p["ncut"] = findall(p["ncut"] .> 0)
  p["cond"] = findall(p["cond"]  .> 0)
  p["spectral"] =findall(p["spectral"]  .> 0)
  p["modularity"] = Int.(p["modularity"]).+1
  # see https://discourse.julialang.org/t/how-to-make-a-named-tuple-from-a-dictionary/10899/14
  return NamedTuple{Tuple(Symbol.(keys(p)))}(values(p))
  #return (;p...) # make  a dictionary into a named tuple
end

function _read_final_with_products(fn::AbstractString)
  G = _read_final(fn)
  pw = _read_products(splitext(fn)[1]*"-products-weighted.json")
  ps = _read_products(splitext(fn)[1]*"-products-simple.json")
  return (G..., products=(simple=ps, weighted=pw), xy=ps.xy)
end


function _read_final_sequence(fn::AbstractString)
  gdata = JSON.parsefile(fn)
  n = gdata["vertices"]
  As = _edgedata_to_sparse.(gdata["sequence"], n)
  # find all instances of edgedata and convert to adjaency matrices...
  return (T=As, dates=Date.(gdata["dates"]), names=string.(gdata["labels"]), orgs=Int.(gdata["orgs"]))
end
