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
## Setup PMC
#=
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
=#
## We can actually verify the maxclique without PMC in this case.
# There's a clique with the same size as the largest core number.
# This is a certificate that you have the max-clique.
# THIS DOES NOT WORK IN GENERAL TO FIND THE MAX-CLIQUE. But since
# we happen to know how to characterize this one, we can prove
# that we have the optimal with some simple code based on the
# maxcore.
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
C1 = simple_clique_heuristic(Tc1)
C2 = simple_clique_heuristic(Tc2)

##
C1 == C2

##
Tcc_ids = C1
##
# This is the largest temporal strong component on the expanded version from Pmc
# Tcc_ids = vec([ 672 728 729 730 733 734 826 735 737 751 792 818 133 80 81 83 88 102 123 124 125 79 137 140 159 161 162 220 229 34 2 3 4 5 7 22 30 32 239 35 36 48 54 67 68 77 536 551 556 559 576 584 508 613 635 639 666 1 677 694 350 241 245 276 278 292 293 294 299 368 375 399 412 480 498 503 795])
## It's actually also a Tcc on the more restrictive version
@assert Tc1[Tcc_ids,Tcc_ids] == ones(length(Tcc_ids),length(Tcc_ids))
##
println.(R1.names[Tcc_ids]);
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
pos2xy(pos::Vector) = [first.(pos) last.(pos)]


##
F = sum(last.(T.T)) # adjacency matrix sum
plot(draw_graph_lines(F,pos2xy(spring(F;iterations=20)))...,
  framestyle=:none, legend=false)

##
#anim = @animate for pos in LayoutIterator(SFDP(;K=0.1,C=0.01,iterations=200), F)
anim = @animate for pos in LayoutIterator(Spring(;initialtemp=0.2,C=0.2,iterations=1000), 10*F .- ones(size(F)...))
  xy = pos2xy(pos)
  plot(draw_graph_lines(F[2:end,2:end],xy[2:end,:])..., marker=:dot, markersize=5, alpha=0.5,
    framestyle=:none, legend=false)
end every 10
gif(anim, "anim.gif", fps=15)

##
function draw_graph_lines_tuple(A::SparseMatrixCSC, xy; shorten::Float64=0)
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
##
plot(draw_graph_lines_tuple(F[2:end,2:end],xy[2:end,:])..., marker=:dot, markersize=5, alpha=0.5,
  framestyle=:none, legend=false)
## Try and do a dynamic layout
using GeometryBasics: Point
using LinearAlgebra
function apply_forces(adj_matrix, old_pos = Vector{Point};
    C::Float64 = 2.0, R::Float64=2.0, temp::Float64=0.025, normalize::Bool=false)

  # The optimal distance bewteen vertices
  N = size(adj_matrix, 1)
  force = similar(old_pos)
  Ftype = eltype(force)
  K = C * sqrt(4.0 / N)
  R = R

  locs = copy(old_pos)
  # Calculate forces
  for i in 1:N
      force_vec = zero(Ftype)
      for j in 1:N
          i == j && continue
          d = norm(locs[j] .- locs[i])
          if adj_matrix[i, j] != zero(eltype(adj_matrix)) || adj_matrix[j, i] != zero(eltype(adj_matrix))
              # F = d^2 / K - K^2 / d
              #F_d = max(adj_matrix[i, j], adj_matrix[j, i])*d / K - R*K^2 / d^2
              # F_d = d / K - K^3 / d^3
              F_d = min(1,max(adj_matrix[i, j], adj_matrix[j, i]))*d / K - R*K^2 / d^2
          else
              # Just repulsive
              # F = -K^2 / d^
              F_d = -R*K^2 / d^2
          end
          # d  /          sin θ = d_y/d = fy/F
          # F /| dy fy    -> fy = F*d_y/d
          #  / |          cos θ = d_x/d = fx/F
          # /---          -> fx = F*d_x/d
          # dx fx
          force_vec += Ftype(F_d .* (locs[j] .- locs[i]))
      end
      force[i] = force_vec
  end
  # Cool down
  #temp = algo.initialtemp / iteration
  # Now apply them, but limit to temperature
  oldnorm = sum(norm.(locs))
  for i in 1:N
      force_mag = norm(force[i])
      scale = min(force_mag, temp) ./ force_mag
      locs[i] += force[i] .* scale
  end
  if normalize
    newnorm = sum(norm.(locs))
    factor = oldnorm/newnorm
    for i in 1:N # rescale to same energy...
      locs[i] *= factor
    end
  end
  return locs
end
function dynamic_layout(T::NamedTuple;
  niter=20,
  niterease=5,
  start=5,
  gamma=0.5,
  droptol=0.01,
  kwargs...
  )
  output = Tuple{Matrix{Float64},SparseMatrixCSC{Float64,Int},Date}[]
  N = size(T.T[1][2],1)
  # start with a circular layout
  pos = [Point{2,Float64}(cos(2*i*pi/N), sin(2*i*pi/N)) for i=1:N]
  Finit = sum(last.(T.T[1:start])) |> x->max.(x,x') # get the initial matrix... make symmetric
  F = min.(Finit,1)
  for init_iter = 1:niter
    pos = apply_forces(F, pos; kwargs...)
  end
  for (date,A) in T.T[start+1:end]
    Fold = copy(F)
    A = max.(A,A')
    A = min.(A,1)
    F = gamma*F + A # update matrix...
    droptol!(F, droptol)
    for ease_iter = 1:niterease
      # interpolate between Fold, F using sin easing
      # From Animations.jl
      # f_sineio(fraction) = sin(pi * fraction - 0.5pi) * 0.5 + 0.5
      amount = sin(pi*ease_iter/(niterease+1)-0.5pi)*0.5 + 0.5
      Fcur = amount*F + (1-amount)*Fold
      droptol!(Fcur, droptol)
      pos = apply_forces(Fcur, pos; kwargs...)
      push!(output, (pos2xy(pos),Fcur,date))
    end
    for iter = 1:niter
      pos = apply_forces(F, pos; kwargs...)
      push!(output, (pos2xy(pos),F,date))
    end
  end
  return output
end
rval = dynamic_layout(T;gamma=0.5,R=1.25,C=1.0,normalize=true)
anim = @animate for (pos,mat,date) in rval
  #@show size(pos)
  matdraw = triu(mat,1)
  plot(draw_graph_lines_tuple(matdraw,pos;shorten=0.1)...,
    line_z = -nonzeros(matdraw)', linewidth=2*nonzeros(matdraw)', alpha=0.5,
    framestyle=:none, legend=false,colorbar=false, clims=(-1,0), size=(800,800))
  scatter!(pos[:,1],pos[:,2], markersize=20, color=T.orgs, alpha=0.5,
    markerstrokewidth=0)
  for i=1:length(T.names)
    annotate!(pos[i,1],pos[i,2], (split(T.names[i], ",")[1], 7, :black))
  end
  title!(string(date))
end every 2
gif(anim, "anim.gif", fps=30)



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
@time cc = HyperModularity.LambdaLouvain(M, zeros(size(M,1)), 0.0)
##
maximum(cc[:,end])
##
function compute_modularity(A::SparseMatrixCSC,c)
    obj = 0
    m = sum(nonzeros(A)/2)
    n = size(A,1)
    # index c...
    C = sparse(1:n, c, 1, n, maximum(c))
    Crv = rowvals(C)
    rv = rowvals(A)
    nz = nonzeros(A)
    for ci=1:size(C,2)
      Cset = Set(Crv[nzrange(C,ci)]) # index the cluster
      for nzci in nzrange(C,ci)
        i = Crv[nzci]
        for nzi in nzrange(A,i)
          if rv[nzi] in Cset
            obj += nz[nzi]
          end
        end
      end
    end
    return obj
end
compute_modularity(M,cc[:,end])
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
