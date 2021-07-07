## Test HyperMoularity.LambdaLouvain with implicit vs. explicit graphs...
# 
include("methods.jl")
G = _build_email_tofrom_graph(data;keepfauci=false,mindegree=2)
##
A = G.A
##
using HyperModularity
cc = HyperModularity.VanillaModularity(Float64.(A),1.0)
##
# k = vec(sum(G.A,dims=2))
# B = sparse(A - (2/sum(A))*(k*k'))
# ccs = HyperModularity.LambdaLouvain(Float64.(B),zeros(size(B,1)),0.0)
##
ccs[:,end]
##
function HyperModularity.ConstructAdj(C::SparseArrays.SparseMatrixCSC,n::Int64)
    """
    ConstructAdj: Construct Adjacency List
    This takes in a sparse adjacency matrix for a graph, and returns an adjacency
    list. While it seems inefficient to store the matrix in multiple ways, as long
    as there are no terrible memory issues, it is VERY handy to be able to quickly
    get a list of neighbors for each node.

    The function also returns the degree vector d. This is NOT the weighted degree,
        d[i] = total number of neighbors of node i
    """

    @show "minenew!"
    rv = C.rowval
    cp = C.colptr
    nz = C.nzval
    Neighbs = Vector{Vector{Int64}}()
    d = zeros(Int64,n)
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        neighs = eltype(rv)[]
        for nzi=cp[i]:cp[i+1]-1
            if nz[nzi] > 0
                push!(neighs, rv[nzi])
            end
        end
        push!(Neighbs,neighs)
        d[i] = length(neighs)
    end

    # d is the number of neighbors. This is the unweighted degree,
    # but note importantly that if the original graph is weighted this is
    # not the same as the degree vector d we will sometimes use
    return Neighbs, d
end

##
k = vec(sum(G.A,dims=2))
B = sparse(A - (1/sum(A))*(k*k'))
ccs = HyperModularity.LambdaLouvain(Float64.(B),zeros(size(B,1)),0.0)
