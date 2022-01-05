##
module PageRank

using SparseArrays, MatrixNetworks, DataStructures

"""
- `maxresidvol::Int` - the maximum residual volume considered, if this is negative,
then we treat it as infinite.

Returns
-------
(x::Dict{Int,Float64},r::Dict{Int,Float64},flag::Int)
"""
function weighted_ppr_push(A::SparseMatrixCSC{T,Int}, seedset,
    alpha::Float64, eps::Float64, maxpush::Int, dvec::Vector{T}, maxresidvol::Int) where T

    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    n = size(A,1)

    x = Dict{Int,Float64}()     # Store x, r as dictionaries
    r = Dict{Int,Float64}()     # initialize residual
    #Q = Int[]        # initialize queue
    Q = Queue{Int}()
    npush = 0.

    if maxresidvol <= 0
        maxresidvol = typemax(Int)
    end

    rvol = 0

    svol = 0.0
    for s in seedset
        svol += dvec[s]
    end
    for s in seedset
        r[s] = dvec[s]/svol
        enqueue!(Q,s)
    end

    pushcount = 0
    pushvol = 0

    @inbounds while length(Q) > 0 && pushcount <= maxpush
        pushcount += 1
        u = dequeue!(Q)

        du = dvec[u] # get the degree

        pushval = r[u] - 0.5*eps*du
        x[u] = get(x,u,0.0) + (1-alpha)*pushval
        r[u] = 0.5*eps*du

        pushval = pushval*alpha

        for nzi in colptr[u]:(colptr[u+1] - 1)
            pushvol += 1
            v = rowval[nzi]
            dv = dvec[v] # degree of v

            rvold = get(r,v,0.)
            if rvold == 0.
                rvol += dv
            end
            rvnew = rvold + pushval*nzval[nzi]/du

            r[v] = rvnew
            if rvnew > eps*dv && rvold <= eps*dv
                #push!(Q,v)
                enqueue!(Q,v)
            end
        end

        if rvol >= maxresidvol
            return x, r, -2
        end
    end

    if pushcount > maxpush
        return x, r, -1, pushcount
    else
        return x, r, 0, pushcount
    end
end

function weighted_ppr_push_solution(A::SparseMatrixCSC{T,Int}, alpha::Float64, seedset, eps::Float64) where T
    maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    dvec = sum(A,dims=2)
    return weighted_ppr_push(A,seed,alpha,eps,maxpush,vec(dvec),0)[1]
end

function weighted_local_sweep_cut(A::SparseMatrixCSC{T,Int}, x::Dict{Int,V}, dvec::Vector{Float64}, Gvol::Float64) where {T,V}
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    n = size(A,1)

    sx = sort(collect(x), by=x->x[2], rev=true)
    S = Set{Int64}()
    volS = 0.
    cutS = 0.
    bestcond = 1.
    beststats = (1,1,1,Gvol-1)
    bestset = Set{Int64}()
    for p in sx
        if length(S) == n-1
            break
        end
        u = p[1] # get the vertex
        #volS += colptr[u+1] - colptr[u]
        volS += dvec[u]

        for nzi in colptr[u]:(colptr[u+1] - 1)
            v = rowval[nzi]
            ew = nzval[nzi]

            if v in S
                cutS -= ew
            else
                cutS += ew
            end
        end
        push!(S,u)
        if cutS/min(volS,Gvol-volS) <= bestcond
            bestcond = cutS/min(volS,Gvol-volS)
            bestset = Set(S)
            beststats = (cutS,min(volS,Gvol-volS),volS,Gvol-volS)
        end
    end
    return bestset, bestcond, beststats
end

function weighted_degree_normalized_sweep_cut!(A::SparseMatrixCSC{T,Int}, x::Dict{Int,V}, dvec::Array{Int}, Gvol::Int) where {T,V}
    colptr = A.colptr
    rowval = A.rowval

    for u in keys(x)
        x[u] = x[u]/dvec[u]
    end

    return weighted_local_sweep_cut(A,x,dvec,Gvol)
end

function weighted_ppr_grow_one(A::SparseMatrixCSC{T,Int},
        seed::Int, alpha::Float64, eps::Float64) where T
    maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    dvec = vec(sum(A,dims=2))
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)
    ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec, 0)[1]
    return weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)
end

"""
acl_diffusion -
G - GraphAndDegrees from SLQ
gamma - sparse matrix regularization
kappa - sparsity regularization

This turns these into alpha, and eps for acl

# This follows from the 1-norm
"""
function acl_diffusion(A,deg,seedset,gamma::T, kappa::T;
        maxpush=nothing, maxresidvol=nothing)  where T
    alpha = 1/(1+gamma)
    #eps = kappa*gamma/alpha
    eps = kappa/(alpha*sum(@view deg[seedset]))
    #
    maxpushval = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    if maxpush != nothing
        maxpushval = round(Int,maxpush)
    end

    ppr = weighted_ppr_push(A,seedset,alpha,eps,maxpushval,deg,
        maxresidvol == nothing ? 0 : maxresidvol)[1]
    x = zeros(size(A,1))
    x[collect(keys(ppr))] .= values(ppr) # almost efficient...
    return x
end

function round_to_cluster(G,x)
    A = G.A
    nnz_dict = Dict{Int,Float64}()
    for (i,v) in enumerate(x)
        if v > 0
            nnz_dict[i] = v
        end
    end
    bestset, bestcond, beststats = weighted_local_sweep_cut(A, nnz_dict, G.deg, sum(G.deg))
    return bestset,bestcond
end

# function degnorm_acl_diffusion(G,seedset, gamma::T, kappa::T;
#     kwargs...) where T
#     x = acl_diffusion(G,seedset,gamma,kappa)
#     x ./= G.deg
#     return x
# end

"""
This computes a full PageRank diffusion with no sparsity
"""
function pr_diffusion(G,seedset,gamma::T)  where T
    alpha = 1/(1+gamma)
    v = sparsevec(seedset, Float64.(G.deg[seedset]), size(G.A,1))
    pr = seeded_pagerank(G.A, alpha, v)
    return pr
end

end
