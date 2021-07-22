using LinearAlgebra
using SparseArrays
using Random


# Expand a hypergraph using the thresholded linear splitting function.
#
#   H = |E| x |V| binary incidence matrix for the hypergraph
#   delta = TL splitting function parameter
# easiest case: delta = 1 is the standard hypergraph cut penalty
#
# order[j] = number of nodes in the j-th hyperedge
#
# This outputs an adjacency matrix A for a directed graph.
#   If H is m x n, then the first n nodes of A correspond to nodes in the original hypergraph
#   the remaining nodes are "auxiliary nodes"
function tl_expansion_inc(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}; delta::Float64=1.0)

        n = size(H,2)
        BigEdges = length(findall(x->x>1,order))
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
            # if order[ee] != nv
            #     @show ee, nv, order[ee], edge
            # end
            if nv == 1
                # ignore
                # println("This")
            # elseif nv == 2
            #     i = edge[1]; j = edge[2]
            #     #A[i,j] += 1; A[j,i] += 1
            #     push!(U,i); push!(V,j); push!(vals,1)
            #     push!(U,j); push!(V,i); push!(vals,1)
            # elseif nv == 3
            #     i = edge[1]; j = edge[2]; k = edge[3]
            #     # A[i,j] += 1/2; A[j,i] += 1/2
            #     # A[k,j] += 1/2; A[j,k] += 1/2
            #     # A[k,i] += 1/2; A[i,k] += 1/2
            #     push!(U,i); push!(V,j); push!(vals,1/2)
            #     push!(U,j); push!(V,i); push!(vals,1/2)
            #     push!(U,i); push!(V,k); push!(vals,1/2)
            #     push!(U,k); push!(V,i); push!(vals,1/2)
            #     push!(U,j); push!(V,k); push!(vals,1/2)
            #     push!(U,k); push!(V,j); push!(vals,1/2)
            else
                # We need to add auxiliary vertices
                for i = edge
                    # A[i,auxpointer] = 1
                    # A[auxpointer+1,i] = 1
                    # A[auxpointer,auxpointer+1] = delta
                    push!(U,i); push!(V,ap); push!(vals,1)
                    push!(U,ap+1); push!(V,i); push!(vals,1)
                end
                push!(U,ap); push!(V,ap+1); push!(vals,delta)
                ap += 2
            end

        end
        A = sparse(V,U,vals,N,N)
        return A
end


function build_augmented_matrix(A,H,gamma,S)
    N = size(H,2)
    d = vec(sum(A,dims=2))
    ei,ej,w = findnz(A)
    nedges = length(ei)
    nedges += 1
    for i = 1:length(S)
        push!(ei,S[i])
        push!(ej,size(A,1)+1)
        push!(w,gamma*d[S[i]])
    end
    Sbar = setdiff(1:N,S)
    for i = 1:length(Sbar)
        push!(ej,Sbar[i])
        push!(ei,size(A,1)+2)
        push!(w,gamma*d[Sbar[i]])
    end
    deg = zeros(size(A,1)+2)
    deg[1:N] = d[1:N]
    return sparse(ei,ej,w,size(A,1)+2,size(A,1)+2),deg
end


using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra

# Controlled growth from a seed set R in a hypergraph. Look at the one hop
# neighborhood, and order all of those nodes by how many hyperedges they are in
# that have nodes from R. Order that last, and take the top k
function TopNeighbors(H::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},R1hop::Vector{Int64},k::Int64)

    if length(R1hop) > k
        # Get all edges touching R
        HR = H[:,R]
        rp = HR.rowval
        edges = unique(rp)

        # Consider how many touch the 1-hop neighborhood
        HL = H[edges,R1hop]

        # For each node in R1hop, compute the number of edges it has that touch R
        d2R = vec(sum(HL,dims=1))

        # order = sortperm(d2R, rev=true)
        b = partialsortperm(d2R, 1:k, rev=true)
        Rmore = R1hop[b]
    else
        Rmore = R1hop
    end

    return union(R, Rmore)
end


# Controlled growth from a seed set R in a hypergraph. Look at the one hop
# neighborhood, and order all of those nodes by what percent of their
# edges touch R. Order that, and take the top k
function BestNeighbors(H::SparseMatrixCSC{Float64,Int64},d::Vector{Float64},R::Vector{Int64},R1hop::Vector{Int64},k::Int64)

    if length(R1hop) > k
        # Get all edges touching R
        HR = H[:,R]
        rp = HR.rowval
        edges = unique(rp)

        # Consider how many touch the 1-hop neighborhood
        HL = H[edges,R1hop]

        # For each node in R1hop, compute the number of edges it has that touch R
        d1 = d[R1hop]
        d2 = vec(sum(HL,dims=1))

        # order = sortperm(d2R, rev=true)
        b = partialsortperm(d2./d1, 1:k, rev=true)
        Rmore = R1hop[b]
    else
        Rmore = R1hop
    end
    return union(R, Rmore)
end


# Simple function for returning ALL the indices where we find a maximum
function findallmax(v)

    l = length(v)
    m = minimum(v)
    M = maximum(v)
    if M == m
        return collect(1:l)
    else
        Inds = Vector{Int64}()
        a,b = findmax(v)
        while a == M
            push!(Inds,b)
            v[b] = m
            a,b = findmax(v)
        end
        return Inds
    end

end


## Delta-Linear (tl = thresholded linear) conductance computation.
# e.g. tl_cond(H,S,d,delta,volA,order)
function tl_cond(H::SparseMatrixCSC,S::Vector{Int64},d::Vector{Float64},delta::Float64,volA::Float64,order::Vector{Int64})

    if volA == 0.0
        volA = sum(d)
    end
    n = length(d)
    volS = sum(d[S])
    cut = tl_cut(H,S,delta,order)

    cond = cut/min(volS, volA-volS)

    return cond, volS, cut

end

## Delta-Linear (thresholded linear) normalized Cut computation.
# e.g. tl_ncut(H,S,d,delta,volA,order)
function tl_ncut(H::SparseMatrixCSC,S::Vector{Int64},d::Vector{Float64},delta::Float64,volA::Float64,order::Vector{Int64})

    if volA == 0.0
        volA = sum(d)
    end
    n = length(d)
    volS = sum(d[S])
    cut = tl_cut(H,S,delta,order)

    cond = cut/min(volS, volA-volS)
    ncut = cut/(volS) + cut/(volA-volS)

    # rncut = round(Int64,ncut)
    # rcut = round(Int64,cut)
    # rcond = round(cond,digits = 4)
    # rvol = round(Int64,volS)

    return cond, ncut, volS, cut

end

# Thresholded linear cut value for a set
# calling e.g. tl_cut(H,S,delta,order)
function tl_cut(H::SparseMatrixCSC{Float64,Int64}, S::Vector{Int64}, delta::Float64,order::Vector{Int64})

    # Check the cut
    HS = H[:,S]
    sumHS = sum(HS,dims = 2)  # Count number of S nodes in each hyperedge
    inds = findall(x->x>0,sumHS)    # Extract hyperedges with > 0 nodes from S
    ES = sumHS[inds]
    verts = order[inds]               # Get the size of these hyperedges

    # Find number of nodes on small side of cut
    SmallSide = round.(Int64,min.(ES, verts-ES))
    # Compute the cardinality-based cut score
    cutval = 0.0
    for j = 1:length(SmallSide)
        sm = SmallSide[j]
        if sm > 0
            if sm < delta
                cutval += sm
            else
                cutval += delta
            end
        end
    end

    return cutval
end


# For a set S in a hypergraph, return the hypergraph local conductance
# score with thresholded linear splitting penalty
# e.g. hlc_tl(H,order,R,S,d,volA,epsilon,delta)
function hlc_tl(H::SparseMatrixCSC{Float64,Int64},order::Vector{Int64},R::Vector{Int64},
    S::Vector{Int64},d::Vector{Float64},volA::Float64,epsilon::Float64,
    delta::Float64,)

    volS = sum(d[S])
    RnS = intersect(R,S)
    volRnS = sum(d[RnS])
    cut = tl_cut(H,S,delta,order)

    lcond = cut/((1+epsilon)*volRnS - epsilon*volS)

    return lcond

end

# Expand a hypergraph using the thresholded linear splitting function.
#
#   Hyperedges = Hyperedge list
#   delta = TL splitting function parameter
function tl_expansion(Hyperedges::Vector{Vector{Int64}}, order::Vector{Int64}, delta::Float64,n::Int64)

        BigEdges = length(findall(x->x>3,order))
        N = n + 2*BigEdges

        ## Build the adjacency matrix
        ap = n+1   # "auxiliary node pointer", points to next "free" aux node

        # Build the sparse matrix
        U = Vector{Int64}()
        V = Vector{Int64}()
        vals = Vector{Float64}()

        for edge = Hyperedges
            nv = length(edge)
            if nv == 2
                i = edge[1]; j = edge[2]
                #A[i,j] += 1; A[j,i] += 1
                push!(U,i); push!(V,j); push!(vals,1)
                push!(U,j); push!(V,i); push!(vals,1)
            elseif nv == 3
                i = edge[1]; j = edge[2]; k = edge[3]
                # A[i,j] += 1/2; A[j,i] += 1/2
                # A[k,j] += 1/2; A[j,k] += 1/2
                # A[k,i] += 1/2; A[i,k] += 1/2
                push!(U,i); push!(V,j); push!(vals,1/2)
                push!(U,j); push!(V,i); push!(vals,1/2)
                push!(U,i); push!(V,k); push!(vals,1/2)
                push!(U,k); push!(V,i); push!(vals,1/2)
                push!(U,j); push!(V,k); push!(vals,1/2)
                push!(U,k); push!(V,j); push!(vals,1/2)
            else
                # We need to add auxiliary vertices
                for i = edge
                    # A[i,auxpointer] = 1
                    # A[auxpointer+1,i] = 1
                    # A[auxpointer,auxpointer+1] = w2
                    push!(U,i); push!(V,ap); push!(vals,1)
                    push!(U,ap+1); push!(V,i); push!(vals,1)
                end
                push!(U,ap); push!(V,ap+1); push!(vals,delta)
                ap += 2
            end

        end
        @show maximum(U), maximum(V), N
        A = sparse(U,V,vals,N,N)
        return A
end

# Given an incidence matrix for a hypergraph and its transpose (having both
# handy makes different parts of the code faster), and a set of nodes R,
# return the immediate neighbors of R that don't include R itself
function get_immediate_neighbors(H::SparseMatrixCSC{Float64,Int64},
    Ht::SparseMatrixCSC{Float64,Int64},R::Vector{Int64})

    Hr = H[:,R]
    rp_r = Hr.rowval
    R_edges = unique(rp_r)

    He = Ht[:,R_edges]
    rp_e = He.rowval
    Rneighbs = unique(rp_e)
    Rn = setdiff(Rneighbs,R)

    return Rn

end

function neighborhood(H::SparseMatrixCSC{Float64,Int64},
    Ht::SparseMatrixCSC{Float64,Int64},R::Vector{Int64})
    Hr = H[:,R]
    rp_r = Hr.rowval
    R_edges = unique(rp_r)

    He = Ht[:,R_edges]
    rp_e = He.rowval
    Rn = unique(rp_e)

    return Rn
end

function neighborlist(H::SparseMatrixCSC{Float64,Int64},
    Ht::SparseMatrixCSC{Float64,Int64})

    Neighbs = Dict()
    n = size(H,2)
    t1 = 0
    t2 = 0
    t3 = 0
    for i = 1:n
        # s = time()
        ivec = H[:,i]
        #n_edges = findnz(ivec)[1]
        n_edges = ivec.nzind    # get neighboring edges
        # t1 += time()-s

        # s = time()
        He = Ht[:,n_edges]      # nodes touching those edges
        rp_e = He.rowval
        neighbs_i = unique(rp_e)
        # t2 += time()-s

        # s = time()
        push!(Neighbs,neighbs_i)
        # t3 += time()-s
    end
    # @show t1, t2, t3
    return Neighbs
end

# Expand a hypergraph using the thresholded linear splitting function.
#
#   H = |E| x |V| binary incidence matrix for the hypergraph
#   delta = TL splitting function parameter
function tl_expansion_inc(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, delta::Float64)

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
            # if order[ee] != nv
            #     @show ee, nv, order[ee], edge
            # end
            if nv == 1
                # ignore
                # println("This")
            elseif nv == 2
                i = edge[1]; j = edge[2]
                #A[i,j] += 1; A[j,i] += 1
                push!(U,i); push!(V,j); push!(vals,1)
                push!(U,j); push!(V,i); push!(vals,1)
            elseif nv == 3
                i = edge[1]; j = edge[2]; k = edge[3]
                # A[i,j] += 1/2; A[j,i] += 1/2
                # A[k,j] += 1/2; A[j,k] += 1/2
                # A[k,i] += 1/2; A[i,k] += 1/2
                push!(U,i); push!(V,j); push!(vals,1/2)
                push!(U,j); push!(V,i); push!(vals,1/2)
                push!(U,i); push!(V,k); push!(vals,1/2)
                push!(U,k); push!(V,i); push!(vals,1/2)
                push!(U,j); push!(V,k); push!(vals,1/2)
                push!(U,k); push!(V,j); push!(vals,1/2)
            else
                # We need to add auxiliary vertices
                for i = edge
                    # A[i,auxpointer] = 1
                    # A[auxpointer+1,i] = 1
                    # A[auxpointer,auxpointer+1] = delta
                    push!(U,i); push!(V,ap); push!(vals,1)
                    push!(U,ap+1); push!(V,i); push!(vals,1)
                end
                push!(U,ap); push!(V,ap+1); push!(vals,delta)
                ap += 2
            end

        end
        # @show maximum(U), maximum(V), length(U), length(V), N, ap
        A = sparse(U,V,vals,N,N)
        return A
end

# Take a list of hyperedges and turn it into a hyperedge incidence matrix
# H. N is the number of nodes in the hypergraph.
#  H(e,u) = 1  iff node u is in hyperedge e
function elist2incidence(Hyperedges::Vector{Vector{Int64}}, N::Int64)
    U = Vector{Int64}()
    E = Vector{Int64}()
    M = length(Hyperedges)
    for enum = 1:length(Hyperedges)
        e = Hyperedges[enum]
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end

    H = sparse(E,U,ones(length(U)),M,N)
    return H
end


# This computes the precision, recall, and F1 score for a set Returned
# compared against a Target set
function PRF(Target,Returned)

    if length(Returned) == 0
        pr = 0; re = 0; F1 = 0
    else
        TruePos = intersect(Returned,Target)
        pr = length(TruePos)/length(Returned)
        re = length(TruePos)/length(Target)
        F1 = 2*(pr*re)/(pr+re)

        if length(TruePos) == 0
            F1 = 0
        end
    end

    return pr, re, F1

end



## Given a binary incidence matrix H for a hypergraph, find the one-hop
# neighborhood of a set of nodes S
function hyper_neighborhood(H::SparseMatrixCSC{Float64,Int64},S::Vector{Int64})

    A = H'*H
    n = size(A,1)
    for i = 1:n
        A[i,i] = 0
    end
    dropzeros!(A)
    return neighborhood(A,S,1)

end

## Given a binary incidence matrix H for a hypergraph, find the one-hop
# neighborhood of a set of nodes S, when considering only hyperedges with
# a maximum number of M nodes
function hyper_neighborhood(H::SparseMatrixCSC{Float64,Int64},S::Vector{Int64},order::Vector{Int64},M::Int64)

    good = findall(x->x<=M,order)
    H = H[good,:]

    ##
    A = H'*H
    n = size(A,1)
    for i = 1:n
        A[i,i] = 0
    end
    dropzeros!(A)
    return neighborhood(A,S,1)

end

## Simple Clique Expansion
#   A[i,j] = number of hyperedges nodes i and j share
function SimpleCliqueExp(H::SparseMatrixCSC{Float64,Int64})

    A = H'*H
    for i = 1:size(A,1)
        A[i,i] = 0.0
    end
    dropzeros!(A)
    return A
end

## Weighted Clique Expansion
#  When performing the clique expansion, for each hyperedge expanded into a
#   clique, multiply each edge in the expansion by 1/order(e)
function WeightedCliqueExpansion(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64})

    m,n = size(H)
    I = Vector{Int64}()
    J = Vector{Int64}()
    vals = Vector{Float64}()
    Hyperedges = incidence2elist(H)
    for e = 1:m
        Edge = Hyperedges[e]
        Ord = order[e]
        for ii = 1:length(Edge)
            for jj = ii+1:length(Edge)
                i = Edge[ii]
                j = Edge[jj]
                push!(I,i); push!(J,j); push!(vals,1/Ord)
            end
        end
        if mod(e,10000)==0
            println("$e")
        end
    end

    A = sparse(I,J,vals,n,n)
    A = sparse(A+A')
    return A
end

function hyper_sweepcut(H::SparseMatrixCSC,x::Vector{Float64},d::Vector{Float64},
    delta::Float64,volA::Float64,order::Vector{Int64};nseeds::Int64=0)
    if volA == 0.0
        volA = sum(d)
    end
    m,n = size(H)
    sorted_indicies = sortperm(-1*x)
    he_included = zeros(m)
    S = Array{Int64,1}([])
    min_cond = 1.0
    min_thd = 0
    volS = 0
    cut = 0
    for i = 1:n
        if x[sorted_indicies[i]] == 0
            break
        end
        push!(S,sorted_indicies[i])
        volS += d[sorted_indicies[i]]
        cut_change = 0
        for k in H.colptr[sorted_indicies[i]]:(H.colptr[sorted_indicies[i]+1]-1)
            j = H.rowval[k]
            he_old = he_included[j]
            he_included[j] += H.nzval[k]
            smallside_old = round(Int64,min(he_old,order[j]-he_old))
            smallside_new = round(Int64,min(he_included[j],order[j]-he_included[j]))
            if smallside_old > 0 && smallside_new == 0
                if smallside_old < delta
                    cut_change -= smallside_old
                else
                    cut_change -= delta
                end
            elseif smallside_old > 0 && smallside_new > 0
                if smallside_old < delta
                    cut_change -= smallside_old
                else
                    cut_change -= delta
                end
                if smallside_new < delta
                    cut_change += smallside_new
                else
                    cut_change += delta
                end
            elseif smallside_old == 0 && smallside_new > 0
                if smallside_new < delta
                    cut_change += smallside_new
                else
                    cut_change += delta
                end
            end
        end
        cut += cut_change
        cond = cut/min(volS, volA-volS)
        # avoid returning extremely small clusters
        if cond < min_cond && i >= nseeds 
            min_cond = cond
            min_thd = i
        end
    end
    return min_cond,sorted_indicies[1:min_thd]
end

function hypergraph_to_bipartite(G)
    m,n = size(G.H)
    Ga = vcat(hcat(spzeros(n,n),G.Ht),hcat(G.H,spzeros(m,m)))
    deg = vec(sum(Ga,dims=2))
    return Ga,deg
end

function CliqueExpansion(G;weighted::Bool=false,binary::Bool=false,thd=50)
    """
    Weighted clique expansion where a hyperedge e is expanded to a
    weighted clique with each edge having weight 1/(|e| - 1)
    """
    n = length(G.deg)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for edge_id = 1:size(G.H,1)
        edge = G.Ht.rowval[G.Ht.colptr[edge_id]:(G.Ht.colptr[edge_id+1]-1)]
        k = length(edge)
        if k > thd
            continue
        end
        if edge_id % 1000 == 0
            @show edge_id,size(G.H,1)
        end
        for i = 1:(k-1)
            ei = edge[i]
            for j = (i+1):k
                ej = edge[j]
                push!(I,ei)
                push!(J,ej)
                if weighted
                    push!(V,1/k)
                else
                    push!(V,1)
                end
            end
        end
    end
    A = SparseArrays.sparse(I,J,V,n,n)
    for i = 1:n; A[i,i] = 0.0; end
    SparseArrays.dropzeros!(A)
    A = SparseArrays.sparse(A+A')
    if binary
        I,J,V = SparseArrays.findnz(A)
        A = SparseArrays.sparse(I, J, 1, n, n)
    end
    deg = vec(sum(A,dims=1))
    return A,deg
end

using PyCall
@pyimport collections as py_collections

# a utility function to read graphs and clusters
function read_dataset(dataset::String;min_size=100,max_size=10000)
    clusters = Dict()
    if dataset == "tripadvisor"
        M = matread("hypergraphs/tripadvisor_H.mat")
        H = M["H"]
        G = LH.graph(H,1.0)
        volA = sum(G.deg)
        Hall = M["H_all"]
        m,n = size(H)
        labelM = zeros(Int64,n,4)
        labelM[:,1] = M["label_code"]
        labelM[:,2] = M["label_region"]
        labelM[:,3] = M["label_locality"]
        labelM[:,4] = M["label_country"]
        zipcodes = M["codes"]
        countries = M["countries"]
        localities = M["localities"]
        regions = M["regions"]
        weights = M["weights"]
        ## Select one type of cluster
        j = 3
        NodeLabels = labelM[:,j]
        if j == 1
            LabelNames = zipcodes
        elseif j == 2
            LabelNames = regions
        elseif j == 3
            LabelNames = localities
        else
            LabelNames = countries
        end
        labels = unique(NodeLabels)
        ## Cluster statistics
        H = Hall
        d = vec(sum(H,dims=1))
        order = round.(Int64,vec(sum(H,dims=2)))
        volA = sum(d)
        m,n = size(H)
        Ht = sparse(H')
        println("Label \t |T| \t cond \t Cluster Name")
        for i = 1:length(labels)
            label = labels[i]
            T = findall(x->x ==label,NodeLabels)
            nT = length(T)
            condS, volS, cutS = tl_cond(H,T,d,1.0,volA,order)
            c = round(condS, digits = 3)
            if nT > 100 && length((LabelNames[label]))> 0 && c < .3
                println("$label \t\t $nT \t $c \t $(LabelNames[label]) \t ")
                clusters[label] = (label,T,LabelNames[label],condS)
            end
        end
    elseif dataset == "amazon"
        M = matread("hypergraphs/AmazonReview5core_H.mat")
        NodeLabels = vec(M["NodeLabels"])
        LabelNames = M["LabelNames"]
        H = M["H"]
        G = LH.graph(H,1.0)
        volA = sum(G.deg)
        node_labels = py_collections.Counter(NodeLabels)
        node_labels = [(key,node_labels[key]) for key in keys(node_labels)]
        node_labels = sort(filter(x->(x[2]>min_size && x[2]<max_size),node_labels),by=x->x[2])
        labels = [x[1] for x in node_labels]
        for i = 1:length(labels)
            label = labels[i]
            T = findall(x->x ==label,NodeLabels)
            nT =length(T)
            condT, volT, cutT = tl_cond(H,T,G.deg,1.0,volA,G.order)
            println("$label \t $nT \t $condT \t $(LabelNames[label])")
            clusters[label] = (label,T,LabelNames[label],condT)
        end
    elseif dataset == "stackoverflow"
        M = matread("hypergraphs/stackoverflow_answer_H.mat")
        LabelMatrix = M["LabelMatrix"]
        LabelNames = M["LabelNames"]
        MainLabels = M["MainLabels"]
        H = M["H"]
        G = LH.graph(H,1.0)
        volA = sum(G.deg)
        for i = 1:length(MainLabels)
            label = MainLabels[i]
            T = findnz(LabelMatrix[:,label])[1]
            nT =length(T)
            condT, volT, cutT = tl_cond(H,T,G.deg,1.0,volA,G.order)
            if condT < 0.2 && 2000 < nT < 10000
                println("$label \t $nT \t $condT \t $(LabelNames[label])")
                clusters[label] = (label,T,LabelNames[label],condT)
            end
        end
    elseif dataset == "walmart"
        M = matread("hypergraphs/Walmart_H.mat")
        NodeLabels = vec(M["NodeLabels"])
        LabelNames = M["LabelNames"]
        H = M["H"]
        G = LH.graph(H,1.0)
        volA = sum(G.deg)
        node_labels = py_collections.Counter(NodeLabels)
        node_labels = [(key,node_labels[key]) for key in keys(node_labels)]
        node_labels = sort(filter(x->(x[2]>100 && x[2]<10000),node_labels),by=x->x[2])
        labels = [x[1] for x in node_labels]
        for i = 1:length(labels)
            label = labels[i]
            T = findall(x->x ==label,NodeLabels)
            nT =length(T)
            condT, volT, cutT = tl_cond(H,T,G.deg,1.0,volA,G.order)
            if condT < 0.5
                println("$label \t $nT \t $condT \t $(LabelNames[label])")
                clusters[label] = (label,T,LabelNames[label],condT)
            end
        end
    elseif dataset == "mathoverflow"
        M = matread("hypergraphs/mathoverflow_answer_H.mat")
        LabelMatrix = M["LabelMatrix"]
        LabelNames = M["LabelNames"]
        MainLabels = M["MainLabels"]
        H = M["H"]
        G = LH.graph(H,1.0)
        volA = sum(G.deg)
        for i = 1:length(MainLabels)
            label = MainLabels[i]
            T = findnz(LabelMatrix[:,label])[1]
            nT =length(T)
            condT, volT, cutT = tl_cond(H,T,G.deg,1.0,volA,G.order)
            if condT < 0.2
                println("$label \t $nT \t $condT \t $(LabelNames[label])")
                clusters[label] = (label,T,LabelNames[label],condT)
            end
        end
    else
        @error("unable to find dataset")
    end
    return H,clusters
end