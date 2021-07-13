## Goal - extract time-varying graphs
include("methods.jl")

## Analysis 1: Temporal strong component.
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
# The result here is a sequence of adjacency matrices, one per day.
T = build_temporal_graphs(data; subset=Tcc_ids)

## Compare the collapsed temporal reachability graph to the rplied to graph
# We remove Fauci from these comparisons.
TG = (T..., A=sum(last.(T.T)) |> A->max.(A,A')) |>
  x-> subset_and_mindegree_cc_filter(x, 2:size(x.A,1))# adjacency matrix sum
TG = TG |> igraph_layout
drawgraph(TG; markersize=15,alpha=0.5,linewidth=2,shownames=true,shorten=0.2)

## Compare the replied to graph to the temporal strong component without Fauci
RT = _build_email_repliedto_graph(data; keepfauci=false) |> igraph_layout
drawgraph(RT; markersize=15,alpha=0.5,linewidth=2,shownames=true,shorten=0.2)

##
setdiff(TG.names, RT.names) # so all RT names are in TG...
##
setdiff(RT.names, TG.names) # so all RT names are in TG...
## Okay, these are actually fairly different graphs!

## Compute the temporal communicability scores
# Peter Grindrod  1 , Mark C Parsons, Desmond J Higham, Ernesto Estrada
# https://doi.org/10.1103/PhysRevE.83.046120
function temporal_communicability(T::NamedTuple, a::Float64)
  Q = prod(map(t->inv(Matrix(I-a*t[2])), T.T))
  broadcast = vec(sum(Q,dims=2))
  recv = vec(sum(Q,dims=1))
  return (Q, broadcast, receive=recv, names=T.names)
end
tcomm = temporal_communicability(T, 0.2)
##
T.names[sortperm(tcomm.broadcast)[1:10]]
##
T.names[sortperm(tcomm.receive)[1:10]]
## Try and do a dynamic layout
using NetworkLayout
pos2xy(pos::Vector) = [first.(pos) last.(pos)]
#= This code is based on the Spring layout in NetworkLayout.jl,
but modified to make it easier to hack. We can also tweak repulsive forces with R,
there is no 'cooling' because we are doign temporal stuff, and also
we can normalize the layout... =#
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
rval = dynamic_layout(T;gamma=0.66,R=1.2,C=0.75,normalize=true,temp=0.005)
##
anim = @animate for (pos,mat,date) in rval
  #@show size(pos)
  matdraw = triu(mat,1)
  plot(draw_graph_lines_vector(matdraw,pos;shorten=0.1)...,
    line_z = -nonzeros(matdraw)', linewidth=2*nonzeros(matdraw)', alpha=0.5,
    framestyle=:none, legend=false,colorbar=false, clims=(-1,0), size=(800,800))
  scatter!(pos[:,1],pos[:,2], markersize=20, color=T.orgs, alpha=0.5,
    markerstrokewidth=0)
  for i=1:length(T.names)
    annotate!(pos[i,1],pos[i,2], (split(T.names[i], ",")[1], 7, :black))
  end
  title!(string(date))
end
gif(anim, "anim.gif", fps=60)

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
@time cc = HyperModularity.LambdaLouvain(M, zeros(size(M,1)), 0.0)[:,end]
##
maximum(cc)
##
function compute_modularity_direct(A::SparseMatrixCSC,c)
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
compute_modularity_direct(M,cc[:,end])
## Assemble the data into a little heatmap...
C = zeros(Int,maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice]
end
# permute by last slice
p= sortperm(C[:,1])
heatmap(C[p,:],color=reverse(theme_palette(:default)[1:8]))
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
## Try and show on our visualization
# This gets VERY hacky because we need to manually align the times.
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
curslice=1
slicedates = first.(T.T)
anim = @animate for (pos,mat,date) in rval
  #@show size(pos)
  # find slice to get communitiies..
  if slicedates[curslice] != date
    print("Finding new slice for $date, from $curslice ", slicedates[curslice])
    while slicedates[curslice] < date
      curslice+=1
    end
    println("to $curslice ", slicedates[curslice])
  end

  matdraw = triu(mat,1)
  plot(draw_graph_lines_vector(matdraw,pos;shorten=0.1)...,
    line_z = -nonzeros(matdraw)', linewidth=2*nonzeros(matdraw)', alpha=0.5,
    framestyle=:none, legend=false,colorbar=false, clims=(-1,0), size=(800,800))
  scatter!(pos[:,1],pos[:,2], markersize=20, color=Int.(C[:,curslice]), alpha=0.5,
    markerstrokewidth=0)
  for i=1:length(T.names)
    annotate!(pos[i,1],pos[i,2], (split(T.names[i], ",")[1], 7, :black))
  end
  title!(string(date))
end
gif(anim, "anim-mod.gif", fps=60)

## Try and figure out cluster 5...
# This assumes there is some real structure there...
# It's around April 20 with Stover as a key member.
println.(findall(C.==5));
println.(collect(zip(findnz(T.T[85][2])[1:2]...) ) .|> x->T.names[[x[1],x[2]]]);
println.(collect(zip(findnz(T.T[86][2])[1:2]...) ) .|> x->T.names[[x[1],x[2]]]);
println.(collect(zip(findnz(T.T[87][2])[1:2]...) ) .|> x->T.names[[x[1],x[2]]]);
println.(collect(zip(findnz(T.T[88][2])[1:2]...) ) .|> x->T.names[[x[1],x[2]]]);

## Nate suggested longer times
rvallong = dynamic_layout(T;gamma=0.66,R=1.2,C=0.75,normalize=true,temp=0.002,niter=40,niterease=10)

##
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
curslice=1
slicedates = first.(T.T)
ndelay = 30 # needs to be set based on the layout... (neaseiter + half niter)
cursliceupdate = 1 # init here to update the very first time...
anim = @animate for (pos,mat,date) in rval[1:600]
  #@show size(pos)
  # find slice to get communitiies..
  # need to be careful here because we want to delay the update
  if cursliceupdate > 0 # we are delaying
    cursliceupdate -= 1
    if cursliceupdate == 0  # apply the update
      while slicedates[curslice] < date
        curslice+=1
      end
    end
  elseif slicedates[curslice] != date && cursliceupdate == 0
    #=print("Finding new slice for $date, from $curslice ", slicedates[curslice])
    while slicedates[curslice] < date
      curslice+=1
    end
    println("to $curslice ", slicedates[curslice], " in ", ndelay, "frames")
    =#
    cursliceupdate = ndelay
  end

  matdraw = triu(mat,1)
  plot(draw_graph_lines_vector(matdraw,pos;shorten=0.1)...,
    line_z = -nonzeros(matdraw)', linewidth=2*nonzeros(matdraw)', alpha=0.5,
    framestyle=:none, legend=false,colorbar=false, clims=(-1,0), size=(800,800))
  scatter!(pos[:,1],pos[:,2], markersize=20, color=Int.(C[:,curslice]), alpha=0.5,
    markerstrokewidth=0)
  for i=1:length(T.names)
    annotate!(pos[i,1],pos[i,2], (split(T.names[i], ",")[1], 7, :black))
  end
  title!(string(date))
end every 2
gif(anim, "anim-mod-long.gif", fps=60)
