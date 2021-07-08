using NetworkLayout
using Plots
#include("../methods.jl")

using Random
import NetworkLayout: IterativeLayout, LayoutIterator
import GeometryBasics: Point

"""
Add a parameter for repulsive force.
"""
struct RepulsiveSpring{Dim,Ptype} <: IterativeLayout{Dim,Ptype}
    C::Float64
    R::Float64
    iterations::Int
    initialtemp::Float64
    initialpos::Vector{Point{Dim,Ptype}}
    seed::UInt
end

function RepulsiveSpring(; dim=2, Ptype=Float64, C=2.0, R=1.0, iterations=100, initialtemp=2.0, initialpos=Point{dim,Ptype}[],
                seed=1)
    if !isempty(initialpos)
        initialpos = Point.(initialpos)
        Ptype = eltype(eltype(initialpos))
        # TODO fix initial pos if list has points of multiple types
        Ptype == Any && error("Please provide list of Point{N,T} with same T")
        dim = length(eltype(initialpos))
    end
    return RepulsiveSpring{dim,Ptype}(C, R, iterations, initialtemp, initialpos, seed)
end

function Base.iterate(iter::LayoutIterator{<:RepulsiveSpring{Dim,Ptype}}) where {Dim,Ptype}
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    N = size(adj_matrix, 1)
    M = length(algo.initialpos)
    rng = MersenneTwister(algo.seed)
    startpos = Vector{Point{Dim,Ptype}}(undef, N)
    # take the first
    for i in 1:min(N, M)
        startpos[i] = algo.initialpos[i]
    end
    # fill the rest with random points
    for i in (M + 1):N
        startpos[i] = 2 .* rand(rng, Point{Dim,Ptype}) .- 1
    end
    # iteratorstate: #iter nr, old pos
    return (startpos, (1, startpos))
end

function Base.iterate(iter::LayoutIterator{<:RepulsiveSpring}, state)
    algo, adj_matrix = iter.algorithm, iter.adj_matrix
    iteration, old_pos = state
    iteration >= algo.iterations && return nothing

    # The optimal distance bewteen vertices
    N = size(adj_matrix, 1)
    force = similar(old_pos)
    Ftype = eltype(force)
    K = algo.C * sqrt(4.0 / N)
    R = algo.R

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
                F_d = d / K - R*K^2 / d^2
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
    temp = algo.initialtemp
    # Now apply them, but limit to temperature
    for i in 1:N
        force_mag = norm(force[i])
        scale = min(force_mag, temp) ./ force_mag
        locs[i] += force[i] .* scale
    end

    return locs, (iteration + 1, locs)
end

anim = @animate for pos in LayoutIterator(RepulsiveSpring(;initialtemp=0.025,C=2,R=2,iterations=1000), F)
  xy = pos2xy(pos)
  #@show extrema.(eachcol(xy))
  plot(draw_graph_lines(F[2:end,2:end],xy[2:end,:])...,
    line_z = nonzeros(F[2:end,2:end]), marker=:dot, markersize=5, alpha=0.5,
    framestyle=:none, legend=false)
end every 10
gif(anim, "anim.gif", fps=15)

##
xy = igraph_layout(F,"fr")
plot(draw_graph_lines(F[2:end,2:end],xy[2:end,:])..., marker=:dot, markersize=5, alpha=0.5,
  framestyle=:none, legend=false)
##
plot([1,2,Inf,3,4,Inf,5,6,Inf],[4,3,Inf,2,1,Inf,5,6,Inf], line_z = [5, 5, 5, -1, -1, -1, 3,3,3],
    linewidth=[0.4 4 5])

##
function draw_graph_lines_tuple(A::SparseMatrixCSC, xy)
  if issymmetric(A)
      ei,ej = findnz(triu(A,1))[1:2]
  else
      ei,ej = findnz(A)[1:2]
  end
  # find the line segments
  lx = Vector{Tuple{Float64,Float64}}[]
  ly = Vector{Tuple{Float64,Float64}}[]
  lx = zeros(0)
  ly = zeros(0)
  for nz=1:length(ei)
      src = ei[nz]
      dst = ej[nz]
      push!(lx, (xy[src,1], xy[dst,1]))
      push!(ly, (xy[src,2], xy[dst,2])
  end
  return lx, ly
end#
