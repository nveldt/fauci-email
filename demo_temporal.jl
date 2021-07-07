## Goal - extract time-varying graphs
include("methods.jl")

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

function _build_tofrom_edges(emails;
            edges = Tuple{Int,Int}[],
            weights = Float64[],
            keepcc::Bool=false,
            keepfauci::Bool=true,
            maxset::Int = typemax(Int))
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

using Printf
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
R1 = temporal_reachability(data;expandstrongcomponents=false)
R2 = temporal_reachability(data)

## Find temporal components
using DelimitedFiles
Tc1 = min.(R1.R,R1.R')
Tc2 = min.(R2.R,R2.R')
##
writedlm("fauci-reachability.txt", zip(findnz(triu(Tc2,1))[1:2]...))
## Run PMC on this to find the largest clique
# dgleich@recurrent:/p/mnt/scratch/corrclus/pmc$ \
#   ./pmc -f ~/Dropbox/research/2021/06-02-fauci-emails/fauci-email/fauci-reachability.txt -v
pmc_output = """
File Name ------------------------ /home/dgleich/Dropbox/research/2021/06-02-fauci-emails/fauci-email/fauci-reachability.txt
workers: 12
self-loops: 0
Reading time 0.0215292
|V|: 1303
|E|: 18444
p: 0.0217435
d_max: 349
d_avg: 28.3101
explicit reduce is set to 4 seconds
K: 77
k-cores time: 0.000306129, ub: 77
*** [pmc heuristic: thread 1]   current max clique = 77,  time = 0.00805092 sec
[pmc heuristic]	 mc = 77
Heuristic found clique of size 77 in 0.013005 seconds
[pmc: heuristic]  Maximum clique: 672 728 729 730 733 734 826 735 737 751 792 818 133 80 81 83 88 102 123 124 125 79 137 140 159 161 162 220 229 34 2 3 4 5 7 22 30 32 239 35 36 48 54 67 68 77 536 551 556 559 576 584 508 613 635 639 666 1 677 694 350 241 245 276 278 292 293 294 299 368 375 399 412 480 498 503 795
Heuristic found optimal solution.
"""
##
# This is the largest temporal strong component on the expanded version
Tcc_ids = vec([ 672 728 729 730 733 734 826 735 737 751 792 818 133 80 81 83 88 102 123 124 125 79 137 140 159 161 162 220 229 34 2 3 4 5 7 22 30 32 239 35 36 48 54 67 68 77 536 551 556 559 576 584 508 613 635 639 666 1 677 694 350 241 245 276 278 292 293 294 299 368 375 399 412 480 498 503 795])
## It's actually also a Tcc on the more restrictive version
@assert Tc1[Tcc_ids,Tcc_ids] == ones(length(Tcc_ids),length(Tcc_ids))
##
println.(R.names[Tcc_ids]);
## Build the sequence of adjacency matrices on this subset of indices
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
T = build_temporal_graphs(data; subset=Tcc_ids)

## Try some fun layouts and movies...
using NetworkLayout


## Multislice community detection
# http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
#=
Matlab code...
Example 2: Ordered Multislice Network
Define the cell array A of square symmetric NxN matrices of equal size each representing one of the T ordered, undirected network "slices". After setting values for the parameters gamma and omega, the multislice modularity with nearest-slice identity arcs of equal strength omega is calculated by
N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);
=#
function expand_to_slicetime_matrix(T::NamedTuple;
  gamma::Float64=1.0,omega::Float64=0.5)
  dates, mats = first.(T.T), last.(T.T)
  nT = length(dates)
  N = size(mats[1],1)
  @assert(all(N .== size.(mats,1)))
  B = spzeros(nT*N, nT*N) # we are gonna be really naughty with sparse construction
  # TODO, make this much better
  twomu = 0
  slices = UnitRange{Int}[]
  for t=1:nT
    mat = mats[t]
    mat = max.(mat,mat') # make it symmetric...
    k = vec(sum(mat;dims=1)); # column sums...
    twom = sum(k)
    twomu += twom
    slice = (1:N) .+ (t-1)*N
    push!(slices,slice)
    B[slice,slice] = mat-gamma*k*k'/twom
  end
  twomu = twomu+2*omega*N*(nT-1)
  B = B + omega*spdiagm(nT*N, nT*N, -N => ones(N*(nT-1)), N=>ones(N*(nT-1)))
  return B, slices
end
M, slices = expand_to_slicetime_matrix(T)
##
using HyperModularity
cc = HyperModularity.LambdaLouvain(M, zeros(size(M,1)), 0.0)
##
maximum(cc[:,end])
##
function exact_temporal_modularity(B)
    #m = sum(nonzeros(A))/2
    #d = vec(sum(A,dims = 2))
    n = size(A,1)
    lam_mod = 1/(2*m)
    c, D =  LazyExactLambdaCC(A,lam_mod,false)
    modularity_score = compute_modularity(A,d,c)
    return c, modularity_score
end

## Assemble the data into a little heatmap...
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
# permute by last slice
p= sortperm(C[:,end])
heatmap(C[p,:],color=reverse(theme_palette(:default)[1:8]))
ylabel!("Nodes in Temporal Strong Component")
xticks!(1:length(T.T), string.(first.(T.T)), xrotation=90)
## Show by orgs
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
# permute by last slice
p= sortperm(C[:,end])
heatmap(T.orgs[p]*ones(1,100),color=reverse(theme_palette(:default)[1:8]))
ylabel!("Nodes in Temporal Strong Component")
xticks!(1:length(T.T), string.(first.(T.T)), xrotation=90)
##
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
# permute by orgs
p= sortperm(T.orgs)
heatmap(C[p,:],color=reverse(theme_palette(:default)[1:8]))
ylabel!("Nodes in Temporal Strong Component")
xticks!(1:length(T.T), string.(first.(T.T)), xrotation=90)

##
function myfunc(T::NamedTuple{Names,Types}) where {Names,Types, in(Names,:A)}
  @show Names
  @show "A in Names"
end
function myfunc(T::NamedTuple{Names,Types}) where {Names,Types, in(Names,:G)}
  @show Names
  @show "G in Names"
end
