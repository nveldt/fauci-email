##
include("../methods.jl")
include("../methods_temporal.jl")
##
T = _read_final_sequence("fauci-email-bydate-sequence-tofrom.json")

## Show the collapsed graph with Fauci
TG = (T..., A=sum(T.T) |> A->max.(A,A')) # adjacency matrix sum
TG = igraph_layout(TG; random=false)
drawgraph(TG; markersize=15,alpha=0.5,
  linewidth=2,shownames=true, markerstrokecolor=:white,
  markerstrokewidth=3,markerstrokealpha=1,shorten=0.1)

## Show the collapsed graph without Fauci
TG = (T..., A=sum(T.T) |> A->max.(A,A')) |>
  x-> subset_and_mindegree_cc_filter(x, 2:size(x.A,1))# adjacency matrix sum
TG = igraph_layout(TG; random=false)
drawgraph(TG; markersize=15,alpha=0.5,
  linewidth=2,shownames=true, markerstrokecolor=:white,
  markerstrokewidth=3,markerstrokealpha=1,shorten=0.1)

## Show the table for temporal communicatbility
aval = 0.02
tcomm = temporal_communicability(T, aval)
##
centrality_types = ["broadcast $(aval)" => "broadcast \$$(aval)\$",
  "receive $(aval)" => "receive \$$(aval)\$"]
topk = 10
results = Dict(map( g-> begin
  s = Symbol(split(g)[1])
  x = tcomm[s]
  p = sortperm(x, rev=true)
  g=>(topk = (tcomm.names[p[1:topk]] .=> x[p[1:topk]]), x=x, names=tcomm.names)
end, first.(centrality_types)))

##
_write_score_table(results,centrality_types)
## Remove the nodes that only have a direct link to Fauci with one email...
# and also have degree 1.
# also remove the first few emails
TG = (T..., A=sum(T.T) |> A->max.(A,A')) # adjacency matrix sum
tokeep = setdiff(1:length(T.names), findall(TG.A[1,:] .== 1))
##
T2 = (T=map(A->A[tokeep,tokeep], T.T[3:end]), names=T.names[tokeep], dates=T.dates[3:end], orgs=T.orgs[tokeep])
TG2 = igraph_layout((A=sum(T2.T), T2.names, T2.orgs); random=false)
drawgraph(TG2; markersize=15,alpha=0.5,
  linewidth=2,shownames=true, markerstrokecolor=:white,
  markerstrokewidth=3,markerstrokealpha=1,shorten=0.1)

##
## Now do the temporal modularity
using HyperModularity
M, slices = expand_to_slicetime_matrix(T2)
@time cc = HyperModularity.LambdaLouvain(M, zeros(size(M,1)), 0.0)[:,end]
## Assemble the data into a little heatmap...
# only show the nodes that are in more than one group
#=
pyplot()
using Measures
C = zeros(Float64,maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice]
end

# make a list of emails by time.
E=zeros(Int, size(C)...)
for i=1:size(C,1)
  E[i,:] .= map(A->sum(findnz(max.(A,A'))[1].==i), T2.T)
end
E = sparse(E)
# find the first appearance for each node and filter by NaN before/after (1 week) that
for i=1:size(C,1)
  firsttime = findfirst(map(A->sum(findnz(max.(A,A'))[1].==i) > 0, T2.T))
  lasttime = findlast(map(A->sum(findnz(max.(A,A'))[1].==i) > 0, T2.T))
  C[i,1:firsttime-1] .= NaN
  C[i,lasttime + 7:end] .= NaN
end
ngroups = length.(unique.(eachrow(C)))
# permute by last slice
#p= sortperm(C[:,end])
#p = findall(ngroups .> 1)
p = sortperm(collect(zip(ngroups, C[:,end])))
plot(color=cgrad(:Set3_7,categorical=true))
heatmap!(C[p,:],color=cgrad(:Set3_7,categorical=true),
  size=(600,600),colorbar=false)
#spy!(E[p,:])
# simulate spy with a scatter...
si,sj = findnz(E[p,:])[1:2]
scatter!(sj,si,markersize=1, markerstrokecolor=:black,
  markerstrokewidth=0.25,  color=:white, label="")
ylabel!("Nodes in Temporal Strong Component with degree >= 2")
# find each monday
monday_ind = findall(dayofweek.(T2.dates) .== 1)
# NOTE, -7 because pyplot rotates in place... ugh., terrible hack.
xticks!(monday_ind .- 7, string.(T2.dates[monday_ind]), xtickfontrotation=45,xtickfontsize=8)
vline!(monday_ind,color=:black,alpha=0.25,label="",linewidth=1.5)
hline!([findfirst(p .== 1)+0.5, findfirst(p .== 1)-0.5],color=:black,label="",linewidth=0.5)
yticks!(1:length(T2.names), T2.names[p], ytickfontsize=7)
=#
##
#savefig("figures/temporal-modularity.pdf")
## Try a variation
# This was better!
pyplot()
# permute by last slice
#p= sortperm(C[:,end])
#p = findall(ngroups .> 1)
heatmap(C[p,:],color=cgrad(:Set3_7,categorical=true),
  size=(600,700),colorbar=false, framestyle=:grid)
# simulate spy with a scatter...
si,sj = findnz(E[p,:])[1:2]
scatter!(sj,si,markersize=1, markerstrokecolor=:black,
  markerstrokewidth=0.25,  color=:white, label="")
#title!("Nodes in Temporal Strong Component with degree >= 2")
# find each monday
monday_ind = findall(dayofweek.(T2.dates) .== 1)
# NOTE, -7 because pyplot rotates in place... ugh., terrible hack.
xticks!(monday_ind .-7, string.(T2.dates[monday_ind]), xtickfontrotation=45,xtickfontsize=8)
xlabel!("Mondays")
vline!(monday_ind,color=:grey,label="",linewidth=0.5)
hline!([findfirst(p .== 1)+0.5, findfirst(p .== 1)-0.5],color=:black,label="",linewidth=0.5)
#yticks!(1:length(T2.names), T2.names[p], ytickfontsize=7)
yticks!(1:length(T2.names), ["     " for _ in 1:length(T2.names)] )
for i=1:size(C,1)
  firsttime = findfirst(E[p[i],:] .> 0)
  if p[i] == 1
    annotate!(firsttime, i, Plots.text(T2.names[p[i]]*" ", 10, :black, :right))
  else
    annotate!(firsttime, i, Plots.text(T2.names[p[i]]*" ", 8, :darkgray, :right))
  end
end
events=[Date("2020-02-26") => "Vice Pres",
        Date("2020-02-28") => "First death",
        Date("2020-03-13") => "Emergency",
        Date("2020-03-06") => "Funding"]
for (ed,name) in events
  annotate!(findfirst(T2.dates .== ed),
           length(T2.names)+1,
           Plots.text(" "*name, 10, :bottom, rotation=90))
end
plot!(left_margin=15mm,top_margin=28mm,bottom_margin=10mm)
##
savefig("figures/temporal-modularity.pdf")
##


## Dynamic layout
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
pos2xy(pos::Vector) = [first.(pos) last.(pos)]
function dynamic_layout(T::NamedTuple;
  niter=20,
  niterease=5,
  start=5,
  gamma=0.5,
  droptol=0.01,
  kwargs...
  )
  output = Tuple{Matrix{Float64},SparseMatrixCSC{Float64,Int},Date}[]
  N = size(T.T[1],1)
  # start with a circular layout
  pos = [Point{2,Float64}(cos(2*i*pi/N), sin(2*i*pi/N)) for i=1:N]
  Finit = sum(T.T[1:start]) |> x->max.(x,x') # get the initial matrix... make symmetric
  F = min.(Finit,1)
  for init_iter = 1:niter
    pos = apply_forces(F, pos; kwargs...)
  end
  for j in (start+1):length(T.T)
    date = T.dates[j]
    A = T.T[j]
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
#rval = dynamic_layout(T2;gamma=0.66,R=1.2,C=0.75,normalize=true,temp=0.005)
rval = dynamic_layout(T2;gamma=0.66,R=1.2,C=0.75,normalize=true,temp=0.002,niter=40,niterease=10)

##
gr()
C = zeros(maximum(slices[1]), length(slices))
for (i,slice) in enumerate(slices)
  C[:,i] = cc[slice,end]
end
curslice=1
slicedates = T2.dates
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

  nodecolor = cgrad(:Set3_7, categorical=true)[Int.(C[:,curslice])]

  matdraw = triu(mat,1)
  plot(draw_graph_lines_vector(matdraw,pos;shorten=0.1)...,
    line_z = -nonzeros(matdraw)', linewidth=2*nonzeros(matdraw)', alpha=0.5,
    framestyle=:none, legend=false,colorbar=false, clims=(-1,0), size=(800,800))
  scatter!(pos[:,1],pos[:,2], markersize=20, color=nodecolor, alpha=1,
    markerstrokewidth=0)
  for i=1:length(T2.names)
    annotate!(pos[i,1],pos[i,2], (split(T2.names[i], ",")[1], 7, :black))
  end
  title!("\n"*string(date))
  plot!(dpi=150)
end
##
#gif(anim, "figures/anim-mod.gif", fps=60)
##
mp4(anim, "figures/anim-mod.mp4"; fps=60)
