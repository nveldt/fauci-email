using HyperModularity
using SparseArrays

function _build_email_hypergraph(data;
    maxset::Int=typemax(Int), mindegree::Int=0, keepfauci=true,
    hyperedgeparts=("sender","recipients","cc"))
  hyperedges = Vector{Vector{Int}}()
  emails=data["emails"]
  names=data["names"]
  orgs=data["clusters"]
  idfauci=findfirst(names .== "fauci, anthony")-1 # offset
  countit = 0
  for thread in emails
    for email in thread
      people = Set{Int}()
      for group in hyperedgeparts
        map(p->push!(people, p), email[group])
      end
      if length(people) <= maxset
        if keepfauci == false
            people = setdiff(collect(people),idfauci)
        end
        # discard hyperedges if they involve 1 or 0 nodes
        @assert(length(unique(people)) == length(people))
        if length(people) > 1
          push!(hyperedges, people.+ 1)
        end
      end
    end
  end
  n = length(names)
  m = length(hyperedges)
  H, HypEdges, weights = convert_to_weighted(hyperedges,ones(m),n)

  m,n = size(H)
  println("Raw data has ")
  println("  $n nodes and $m hyperedges")

  # filter out low degree nodes
  d = vec(sum(H;dims=1))
  d = hypergraph_sum_degree(H,weights)
  filtnodes = d .>= mindeg
  Hd = H[:,filtnodes]
  filtids = findall(filtnodes) # get the ids from original dataset
  namesd = names[filtids]
  orgsd = orgs[filtids]

  # combine identical edges and reweight
  hyperedges = HyperModularity.incidence2elist(Hd)
  H2, HypEdges2, weights2 = convert_to_weighted(hyperedges,weights,size(Hd,2))

  # Take largest connected component
  m,n = size(H2)
  A = [SparseArrays.spzeros(n,n) H2'; H2 SparseArrays.spzeros(m,m) ]
  Acc,ccfilt = largest_component(A)
  ccfilt_nodes = ccfilt[1:n]
  ccfilt_edges = ccfilt[n+1:end]

  # nodeids = filtids[ccfilt_nodes]
  names_final = namesd[ccfilt_nodes]
  orgs_final = orgsd[ccfilt_nodes]
  Hfinal = H2[ccfilt_edges,ccfilt_nodes]
  weights_final = weights2[ccfilt_edges]
  edgelist = incidence2elist(Hfinal)

  return (H=Hfinal, elist=edgelist, orgs = orgs_final,names = names_final,weights=weights_final)
end

function unroll_weighted_hypergraph(H::NamedTuple)
  # replace an integer weighted hypergraph with a hypergraph
  # where each hyperedge is listed k times if it has weight k
  # in the original hypergraph.
  # Seems silly, but need this to run certain hypergraph PageRank
  # algorithms...
  n = size(H.H,2)
  Elist = incidence2elist(H.H)
  E2 = Vector{Vector{Int64}}()
  for e = 1:length(Elist)
      num = H.weights[e]
      for k = 1:num
          push!(E2,Elist[e])
      end
  end
  H2 = elist2incidence(E2,n)
  return H2
end

function project_hypergraph(H::NamedTuple; distribute_hyperedge::Bool=false )
  # Clique expansion
  n = length(H.names)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e = 1:length(H.elist)
    edge = H.elist[e]
    if maximum(edge) > n
      @show edge
    end
    clique_edge_scale = 1/binomial(length(edge),2)
    for i = 1:length(edge)
      for j = i+1:length(edge)
        push!(I,edge[i])
        push!(J,edge[j])
        if distribute_hyperedge
          push!(V,H.weights[e]*clique_edge_scale)
        else
          push!(V,H.weights[e])
        end
      end
    end
  end
  A = sparse(I,J,V,n,n)
  A = dropzeros!(max.(A,A'))
  return A
end

function convert_to_weighted(hyperedges,weights,n)
  HypDict = Dict()
  for eid = 1:length(hyperedges)
    e = hyperedges[eid]
    if haskey(HypDict,e)
      HypDict[e] += weights[eid]
    else
      HypDict[e] = weights[eid]
    end
  end

  HypEdges = Vector{Vector{Int64}}()
  new_weights = Float64[]
  for edge in keys(HypDict)
    push!(HypEdges,edge)
    push!(new_weights,HypDict[edge])
  end

  H = HyperModularity.elist2incidence(HypEdges,n)
  return H, HypEdges, new_weights
end

function hypergraph_sum_degree(H,weights)
  # degree that you get from a clique expansion
  n = size(H,2)
  sum_d = zeros(n)
  order = vec(sum(H;dims = 2))
  for i = 1:n
    edge_ids = findall(x->x>0,H[:,i])
    for j in edge_ids
      sum_d[i] += weights[j]*(order[j]-1)
    end
  end
  return sum_d
end

function hyperstcut(H,s,t,weights,delta = 1.0)
  order = vec(round.(Int64,sum(H;dims = 2)))
  n = size(H,2)
  @assert(s <= n)
  @assert(t <= n)
  A = tl_expansion_inc(H,order,weights,delta)
  F = maxflow(Float64.(A),s,t,0.0)
  S = intersect(source_nodes_min(F),collect(1:n))

  return S
end

##
function plot_star_expansion(H)
  m,n = size(H)
  A = [SparseArrays.spzeros(n,n) H'; H SparseArrays.spzeros(m,m)]
  xy = igraph_layout(A)
  p = plot(draw_graph_lines_vector(A,xy;shorten=0.0);label=false, alpha=0.1,framestyle=:none, linecolor=1)
  scatter!(p,H.xy[1:n,1],H.xy[1:n,2])
  scatter!(p,H.xy[n+1:n+m,1],H.xy[n+1:n+m,2],color = :blue, markerstrokewidth = 0.0)
end

# Expand a hypergraph using the thresholded linear splitting function.
#
#   H = |E| x |V| binary incidence matrix for the hypergraph
#   delta = thresholded-linear splitting function parameter
  function tl_expansion_inc(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, weights::Vector{Float64},delta::Float64)

    n = size(H,2)
    BigEdges = length(findall(x->x>3,order))
    N = n + 2*BigEdges

    Hyperedges = incidence2elist(H)

    ## Build the adjacency matrix
    ap = n+1   # "auxiliary node pointer", points to next "free" aux node

    # Build the sparse matrix
    U = Vector{Int64}()
    V = Vector{Int64}()
    vals = Vector{Float64}()

    for ee = 1:length(Hyperedges)
        edge = Hyperedges[ee]
        nv = length(edge)
        wts = weights[ee]
        # if order[ee] != nv
        #     @show ee, nv, order[ee], edge
        # end
        if nv == 1
            # ignore
            # println("This")
        elseif nv == 2
            i = edge[1]; j = edge[2]
            #A[i,j] += 1; A[j,i] += 1

            push!(U,i); push!(V,j); push!(vals,wts)
            push!(U,j); push!(V,i); push!(vals,wts)
        elseif nv == 3
            i = edge[1]; j = edge[2]; k = edge[3]
            # A[i,j] += 1/2; A[j,i] += 1/2
            # A[k,j] += 1/2; A[j,k] += 1/2
            # A[k,i] += 1/2; A[i,k] += 1/2
            push!(U,i); push!(V,j); push!(vals,wts/2)
            push!(U,j); push!(V,i); push!(vals,wts/2)
            push!(U,i); push!(V,k); push!(vals,wts/2)
            push!(U,k); push!(V,i); push!(vals,wts/2)
            push!(U,j); push!(V,k); push!(vals,wts/2)
            push!(U,k); push!(V,j); push!(vals,wts/2)
        else
            # We need to add auxiliary vertices
            for i = edge
                # A[i,auxpointer] = 1
                # A[auxpointer+1,i] = 1
                # A[auxpointer,auxpointer+1] = delta
                push!(U,i); push!(V,ap); push!(vals,wts)
                push!(U,ap+1); push!(V,i); push!(vals,wts)
            end
            push!(U,ap); push!(V,ap+1); push!(vals,delta*wts)
            ap += 2
        end

    end
    # @show maximum(U), maximum(V), length(U), length(V), N, ap
    A = sparse(U,V,vals,N,N)
    return A
end


##
function _read_final_hypergraph(fn::AbstractString)
  hdata = JSON.parsefile(fn)
  n = hdata["vertices"]
  m = hdata["hyperedges"]
  nmaps = hdata["incidences"]
  X = reshape(Int.(hdata["hyperedgedata"]), 2, nmaps)
  H = sparse(Int.(X[1,:]).+1, Int.(X[2,:]).+1, 1.0, m, n)
  weights = Float64.(hdata["weights"])
  return (H=H,
    elist=incidence2elist(H),
    orgs=Int.(hdata["orgs"]),
    names=string.(hdata["labels"]),
    weights=weights,)
end
