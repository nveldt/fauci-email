using Dates
function _simple_check_duplicate(emails, email)
  for e in emails
    if e["sender"] == email["sender"] &&
      e["recipients"] == email["recipients"] &&
      e["cc"] == email["cc"]
      return true
    end
  end
  return false
end
""" get a dictionary that lists emails for each day. """
function _emails_by_day_and_time(data;keepduplicates::Bool=false)
  dt = DateFormat("yyyy-mm-ddTHH:MM:SS")
  emails_by_day = Dict{Date, Vector{Any}}()
  emails_by_time = Dict{DateTime, Vector{Any}}() # use to check for duplicates
  duplicates=0
  for thread in data["emails"]
    for email in thread

      t = DateTime(email["time"][1:19], dt)
      day = Date(t)
      if haskey(emails_by_time,t)
        if keepduplicates==false &&
          _simple_check_duplicate(emails_by_time[t], email)==true
          duplicates += 1
          continue # skip this email... it's a duplicate
        else
          push!(emails_by_time[t], email)
        end
      else
        emails_by_time[t] = [email]
      end

      if haskey(emails_by_day, day)
        push!(emails_by_day[day], email)
      else
        emails_by_day[day] = [email]
      end
    end
  end

  # give back a sorted list of days and emails by day
  byday = sort!(collect(emails_by_day), by=x->x[1])
  bytime = sort!(collect(emails_by_time), by=x->x[1])

  return (bydate=byday, bytime=bytime, names=data["names"], orgs=data["clusters"],
    keepduplicates, duplicates)
end

function _vector_of_sets_to_sparse(adjlist::Vector{Set{Int}};
    n=maximum(Iterators.flatten(adjlist)),
    weights=nothing)
  # just copy out edges
  edges = Tuple{Int,Int}[]
  for (ei,neighs) in enumerate(adjlist)
    for ej in neighs
      push!(edges, (ei,ej))
    end
  end
  return sparse(first.(edges),last.(edges),weights === nothing ? 1 : weights, n, n)
end


""" Build the temporal reachability matrix, this means
that i->j if there is a temporal path from i to j. Temporal paths are resolved
down to the level of days. So if i and j exchange email on the same day, then
that event is not ordered and can we treat it as a non-directional edge. """
function temporal_reachability(data; keepcc::Bool=true,
    keepfauci::Bool=true, expandstrongcomponents::Bool=true)
  tdata = _emails_by_day_and_time(data)
  n = length(tdata.names)
  @assert issorted(Iterators.map(first, tdata.bydate)) # check sorting
  #R = length(tdata) |> x->zeros(x,x) # allocate a matrix of zeros
  R = [Set{Int}((j,)) for j=1:length(tdata.names)] # allocate reachability sets with self-reachability
  for (date,emails) in reverse(tdata.bydate) # need to do these in reverse order of time.
    edges = _build_tofrom_edges(emails;keepcc,keepfauci)[1]
    Rcopy = deepcopy(R) # need to make a copy because these aren't temporally ordered...
    # TODO, I think you can just use two copies, since the data only gets bigger...
    for (i,j) in edges
      # i emails j, so extend the temporal path from i to everyone j can reach...
      for k in R[j]
        push!(Rcopy[i], k)
      end
    end

    if expandstrongcomponents
      scc_info = scomponents(sparse(first.(edges),last.(edges),1,n,n))
      for (ci,componentsize) in enumerate(scc_info.sizes)
        if componentsize > 2 # if component size is <= 2, we already handled it above since edges must exist
          verts = findall(scc_info.map .== ci)
          for (vi,vj) in Iterators.product(verts, verts)
            for k in R[vj]
              push!(Rcopy[vi], k)
            end
          end
        end
      end
    end
    R = Rcopy # swap...
  end

  # let's turn R into a matrix now...
  return (R=_vector_of_sets_to_sparse(R; n=n), tdata...)
end

function simple_clique_heuristic(A::SparseMatrixCSC)
  cns = corenums(A)[1]
  maxcore = maximum(cns)
  maxcoreids = findall(cns .== maxcore)
  C = greedy_clique(A[maxcoreids,maxcoreids])
  if length(C) != maxcore
    @warn("not maximum clique, so we don't have the largest temporal strong component")
  end
  return maxcoreids[collect(C)]
end
function greedy_clique(A::SparseMatrixCSC)
  B = copy(A)
  fill!(B.nzval, 1)
  B = B - Diagonal(B)
  dropzeros!(B)
  @assert issymmetric(B)
  d = vec(sum(B;dims=2))
  # C is the current clique, F is the set of feasible vertices to add
  # we pick the one with the largest degree.
  function _expand_clique(C::Set{Int}, F::Set{Int})
    maxd = 0
    argmaxF = 0
    for v in F
      if d[v] > maxd
        argmaxF = v
        maxd = d[v]
      end
    end
    push!(C,argmaxF)
    # update F...
    # remove everything in F that doesn't have a link to v.
    filter!(u->B[u,argmaxF] == 1, F)
    if length(F) > 0
      return _expand_clique(C, F)
    else
      return C
    end
  end
  # start that on the vertex of max-degree
  C = Set((argmax(d),))
  F = Set(findnz(B[:,first(C)])[1])
  return _expand_clique(C,F)
end

function build_temporal_graphs(data;
    keepcc::Bool=true,keepfauci::Bool=true,
    emailweight::Bool=false,keepduplicates::Bool=false,
    subset=nothing,keepempty=false)
  tdata = _emails_by_day_and_time(data;keepduplicates)
  nfull = length(tdata.names)
  if subset !== nothing
    sids = sort!(unique(subset))
    #idset =  Set(sids)
    n = length(sids)
  else
    n = nfull
    sids = Colon
  end
  SparseType = emailweight ? SparseMatrixCSC{Int,Float64} : SparseMatrixCSC{Int,Int}
  T = Pair{Date,SparseType}[]
  for (date,emails) in tdata.bydate
    edges,weights = _build_tofrom_edges(emails; keepcc, keepfauci)
    # TODO, could make this more efficient, but, ugh...
    A = sparse(first.(edges),last.(edges), emailweight ? weights : 1, nfull, nfull)
    A = dropzeros!(A[sids, sids]) # filter to subset
    if keepempty
      push!(T,date => A)
    elseif nnz(A) > 0
      push!(T,date => A)
    end
  end

  return (T=T, names=tdata.names[sids], orgs=tdata.orgs[sids],
    keepcc, keepfauci, emailweight, keepempty, keepduplicates,
    duplicates=tdata.duplicates)
end
