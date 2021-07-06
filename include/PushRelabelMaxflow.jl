using MatrixNetworks
using SparseArrays

# Push Relabel solver for maximum s-t flow, minimum s-t cut problems

mutable struct stFlow
    flowvalue::Float64 # gives you the max-flow value
    cutvalue::Float64 # gives min-cut value, which should equal flowvalue,
                      # but may differ by a small tolerance value.
    source_nodes::Vector{Int64} # give the indices of the nodes attached to the source
    C::SparseMatrixCSC # gives the original capacity matrix
    F::SparseMatrixCSC # gives the values of the flows on each edge
    s::Int64  # index of source node
    t::Int64 # index of sink node
end

"""
maxflow

Given a sparse matrix A representing a weighted and possibly directed graph,
a source node s, and a sink node t, return the maximum s-t flow.

flowtol = tolerance parameter for whether there is still capacity available on
            an edge. Helps avoid rounding errors. Default is 1e-6.

Returns F, which is of type stFlow.
"""
function maxflow(B::Union{SparseMatrixCSC,MatrixNetwork},s::Int,t::Int, flowtol::Union{Float64,Int}= 1e-6)

    if flowtol >= .1
        println("flowtol is a tolerance parameter for rounding small residual capacity edges to zero, and should be much smaller than $flowtol. Changing it to default value 1e-6")
        flowtol = 1e-6
    end

    # The code actually assumes a SparseMatrixCSC input
    if typeof(B) <: SparseMatrixCSC
    else
        B = sparse(B)
    end

    N = size(B,1)

    # Extract weights from source s to non-terminal nodes,
    # and from non-terminal nodes to sink node t
    sWeights = Array(B[s,:])
    tWeights = Array(B[:,t])
    NonTerminal = setdiff(collect(1:N),[s t])

    sWeights = sWeights[NonTerminal]
    tWeights = tWeights[NonTerminal]

    # Extract the edges between non-terminal nodes
    A = B[NonTerminal,NonTerminal]

    # A = the matrix of capacities for all nodes EXCEPT the source and sink
    # sWeights = a vector of weights for edges from source to non-terminal nodes
    # tWeights = vector of weights from non-terminal nodes to the sink node t.

    # This is the map from the original node indices to the rearranged
    # version in which the source is the first node and the sink is the last
    Map = [s; NonTerminal; t]

    # Directly set up the flow matrix
    C = [spzeros(1,1) sparse(sWeights') spzeros(1,1);
         sparse(sWeights) A sparse(tWeights);
         spzeros(1,1) sparse(tWeights') spzeros(1,1)]

    # Allocate space for the flow we will calculate
    # In a flow problem, we will eventually need to send flow the reverse
    # direction, so it's important to allocate space for F[i,j] if C[j,i] is an
    # edge, even if C[i,j] is not directed
    Cundir = C+C'
    F = SparseMatrixCSC(N,N,Cundir.colptr,Cundir.rowval,zeros(length(Cundir.rowval)))
    ExcessNodes = vec(round.(Int64,findall(x->x!=0,sWeights).+1))

    # Initialize the Preflow and the excess vector
    for v = ExcessNodes
        F[1,v] = C[1,v]
        F[v,1] = -C[1,v]
    end
    excess = [0;sWeights;0]
    source_nodes, FlowMat, value = Main_Push_Relabel(C,F,ExcessNodes,excess,flowtol)

    smap = sortperm(Map)
    F = stFlow(value, value, sort(Map[source_nodes]),C[smap,smap],FlowMat[smap,smap],s,t)
    return F
end

"""
This maxflow code assumes that A represents the adjacencies between
non-terminal nodes. Edges adjecent to source node s and sink node t
are given by vectors svec and tvec.

This code sets s as the first node, and t as the last node.
"""
function maxflow(A::Union{SparseMatrixCSC,MatrixNetwork},svec::Vector{Float64},tvec::Vector{Float64}, flowtol::Union{Float64,Int}= 1e-6)

    if flowtol >= .1
        println("flowtol is a tolerance parameter for rounding small residual capacity edges to zero, and should be much smaller than $flowtol. Changing it to default value 1e-6")
        flowtol = 1e-6
    end
    if typeof(A) <: SparseMatrixCSC
    else
        A = sparse(A)
    end


    # Directly set up the flow matrix
    C = [spzeros(1,1) sparse(svec') spzeros(1,1);
         sparse(svec) A sparse(tvec);
         spzeros(1,1) sparse(tvec') spzeros(1,1)]

    N = size(C,1)

    # Allocate space for the flow we will calculate
    # In a flow problem, we will eventually need to send flow the reverse
    # direction, so it's important to allocate space for F[i,j] if C[j,i] is an
    # edge, even if C[i,j] is not directed.
    Cundir = C+C'
    F = SparseMatrixCSC(N,N,Cundir.colptr,Cundir.rowval,zeros(length(Cundir.rowval)))
    ExcessNodes = vec(round.(Int64,findall(x->x!=0,svec).+1))

    # Initialize the Preflow and the excess vector
    for v = ExcessNodes
        F[1,v] = C[1,v]
        F[v,1] = -C[1,v]
    end
    excess = [0;svec;0]
    source_nodes, FlowMat, value = Main_Push_Relabel(C,F,ExcessNodes,excess,flowtol)

    F = stFlow(value,value,source_nodes,C,FlowMat,1,N)
end

maxflow(A::Union{SparseMatrixCSC,MatrixNetwork},svec::Vector{Int64},tvec::Vector{Int64},flowtol::Union{Float64,Int}= 1e-6) =
    maxflow(A,float(svec),float(tvec),flowtol)


flow(F::stFlow) =
    F.flowvalue

"""
Given a flow, stored in an stFlow object, return the set of nodes attached to
the source
"""
function source_nodes(F::stFlow,flowtol::Union{Float64,Int}= 1e-6)
    # Run a bfs from the sink node. Anything with distance
    # n is disconnected from the sink. Thus it's part of the minimium cut set
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,flowtol, F.t)
    S = Vector{Int64}()
    for i = 1:n
        if finalHeight[i] == n
            push!(S,i)
        end
    end

    # Sanity checks: source node is on source side, sink node is on sink side
    @assert(~in(F.t,S))
    @assert(in(F.s,S))

    return S
end

# Get the smallest source-side set
function source_nodes_min(F::stFlow,flowtol::Union{Float64,Int}= 1e-6)
    # Run a bfs from the source node. Anything with distance
    # <n is connected to the source. Thus it's part of the minimium cut set
    n = size(F.C,2)
    finalHeight = relabeling_bfs(SparseMatrixCSC(F.C'),SparseMatrixCSC(F.F'),flowtol,F.s)
    S = Vector{Int64}()
    for i = 1:n
        if finalHeight[i] < n
            push!(S,i)
        end
    end

    # Sanity checks: source node is on source side, sink node is on sink side
    @assert(~in(F.t,S))
    @assert(in(F.s,S))

    return S
end

"""
Given a flow, stored in an stFlow object, return the set of nodes attached to
the sink
"""
function sink_nodes(F::stFlow,flowtol::Union{Float64,Int}= 1e-6)
    # Run a bfs from the sink node. Anything with distance < n is sink-attached.
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,flowtol,F.t)
    T = Vector{Int64}()
    for i = 2:n
        if finalHeight[i] < n
            push!(T,i)
        end
    end

    # Sanity checks
    @assert(in(F.t,T))
    @assert(~in(F.s,T))

    return T
end

"""
Gives the cut as a list of edges.
"""
function cut_edges(F::stFlow,flowtol::Union{Float64,Int}= 1e-6)
    # Run a bfs from the sink node to get source and sink sets
    n = size(F.C,2)
    finalHeight = relabeling_bfs(F.C,F.F,flowtol,F.t)
    T = Vector{Int64}()
    S = Vector{Int64}()
    for i = 1:n
        if finalHeight[i] < n
            push!(T,i)
        else
            push!(S,i)
        end
    end

    I,J,V = findnz(F.C[S,T])
    return [S[I] T[J]]
end


"""
Gives the non-terminal cut edges.
"""
function cut_edges_nonterminal(F::stFlow,flowtol::Union{Float64,Int}= 1e-6)
    # Run a bfs from the sink node to get source and sink sets
    Edges = cut_edges(F)
    T = Vector{Int64}()
    S = Vector{Int64}()
    for i = 1:size(Edges,1)
        I = Edges[i,1]
        J = Edges[i,1]
        if I != F.t && I!= F.s && J != F.t && J != F.s
            push!(S,I)
            push!(T,J)
        end
    end
    return [S T]
end

# Main_Push_Relabel returns a maximum flow F and the min s-t cut set S for the
# flow graph C.
#
# C = the capacity matrix for the flow problem.
#   This code assumes node 1 is the source, and node n is the sink.
#   the preflow immediately pushes all flow from the source to create an
#   excess on nodes in the graph.
#
# F = an initial flow. It can be initialize to zero.
#
# excess = the vector of excess values at the start of the algorithm. If F = 0,
#   this is the vector of edge capacities from the implicit source to the graph.
#   If F != 0, then it's the excess from a previous run of the algorithm
function Main_Push_Relabel(C::SparseMatrixCSC,
    F::SparseMatrixCSC,ExcessNodes::Array{Int64},excess::Array{Float64},flowtol::Union{Float64,Int}= 1e-6)

    # here, n includes only one terminal node, the sink
    n = size(C,1)

    height = zeros(Int64,n)      # label/height of each node
    inQ = zeros(Bool,n)          # list whether or not nodes are in the queue

    # Store adjacency list. Because flow can be sent either direction on an
    # arc during the course of the algorithm, it's important to list all neighbors
    # or each node, counting both incoming and outgoing edges
    Neighbs,d = ConstructAdj(C+C',n)

    # We will maintain a queue of active nodes.
    #   An actual queue implementation is available in the DataStructures.jl
    #   Julia package. The performane is nearly identical (and in some cases
    #   slightly slower), thus to minimize dependency on outside packages, we
    #   just use a Vector, rather than an actual implemented Queue
    Queue = Vector{Int64}()

    # Start by saturating edges from source to its neighbors
    # All nodes with nonzero excess are the first to be processed
    for v = ExcessNodes
        push!(Queue,v)
    end
    inQ[ExcessNodes] .= true

    # count the number of nodes that have been relabeled
    relabelings::Int64 = 0

    height = relabeling_bfs(C,F,flowtol,n)  # compute initial distance from sink
    # In the code and comments, height = distance from sink = label of node

    # Continue until the queue no longer contains any active nodes.
    while length(Queue) > 0

        u = pop!(Queue)     # Select a new active node

        inQ[u] = false      # Take it out of the queue

        # discharge flow through node u
        relabelings += discharge!(C,F,Queue,u,Neighbs[u],height,excess,n,d[u],inQ,flowtol)

        # if u is still active, put it back into the queue
        if excess[u] > flowtol
            prepend!(Queue,u)
            inQ[u] = true
        end

        # Global relabeling heuristic for push-relabel algorithm.
        # This periodically recomputes distances between nodes and the sink
        if relabelings == n
            relabelings = 0
            dist = relabeling_bfs(C,F,flowtol)
            height = dist
        end

    end

    # Compute final distances from sink using BFS. Anything with distance
    # n is disconnected from the sink. Thus it's part of the minimium cut set
    finalHeight = relabeling_bfs(C,F,flowtol,n)
    S = Vector{Int64}()
    push!(S,1)          # Include the source node
    for i = 2:n
        if finalHeight[i] == n
            push!(S,i)
        end
    end

    mflow = excess[n]     # the excess at the sink equals the maximum flow value

    return S, F, mflow

end

# Discharege operation: pushes flow away from node u across admissible edges.
# If excess[u] > 0 but no admissible edges exist, we relabel u.
function discharge!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    excess::Array{Float64},n::Int64,du::Int64,inQ::Array{Bool},
    flowtol::Union{Float64,Int}= 1e-6)

    vLocal::Int64 = 1           # Start at the first neighbor of node u
    hu = height[u]
    relabeled = 0

    # As long as there is excess at node u and there is another neighbor to explore...
    while excess[u] > flowtol && vLocal <= du

            # ...grab the next neighbor of node u
            v = uNeighbs[vLocal]

            # ... if edge (u,v) is admissible, push more flow.
            # Otherwise, move to the next neighbor of u
            if hu > height[v] && C[u,v] - F[u,v] > flowtol
                pushflow!(C,F,Queue,u,v,excess,height,inQ,n)
                vLocal += 1
            else
                vLocal += 1
            end
    end

    # if we needed to visit every neighbor of u, we must relabel u,
    # so that at least one admissible edge is created
    if vLocal > du
        relabeled = 1
        relabel!(C,F,Queue,u,uNeighbs,height,du,n,flowtol)
    end

    return relabeled
end

# Relabel sets the label/height of node u to be equal to the minimum label
# such that an admissible edge exists. An edge (u,v) is admissible if
# height[u] = height[v] + 1
function relabel!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    du::Int64,n::Int64,flowtol::Union{Float64,Int}= 1e-6)
   # find smallest new height making a push possible, if such a push is possible

   min_height = Inf
   # search through the neighbors of u
   # and relabel so that height[u] = height[v] + 1 for some v in the neighborhood
   for vLocal = 1:du
       v = uNeighbs[vLocal]
       if C[u,v] - F[u,v] > flowtol
           min_height = min(min_height, height[v])
           height[u] = min_height + 1
       end
   end

end

# Push flow from an active node u to a node v via an admissible edge (u,v)
function pushflow!(C::SparseMatrixCSC,F::SparseMatrixCSC,
    Queue::Vector{Int},u::Int64,v::Int64,excess::Array{Float64},height::Array{Int64},
    inQ::Array{Bool},n::Int64)

    send = min(excess[u], C[u,v] - F[u,v])
    F[u,v] += send
    F[v,u] -= send
    excess[u] -= send
    excess[v] += send

    # If v isn't in the queue, isn't the sink, isn't the source,
    # and is active, then add it to the Queue
    if ~inQ[v] && v < n && v > 1
        prepend!(Queue,v)
        inQ[v] = true
    end
end

# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    d = zeros(Int64,n)
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
        d[i] = ci[i+1]-ci[i]
    end

    # d is the number of neighbors. This is the unweighted degree,
    # but note importantly that if the original graph is weighted this is
    # not the same as the degree vector d we will sometimes use
    return Neighbs, d

end

# Given initial capacity matrix C and flow matrix F, compute the distance
# from each node to the specified "start" node.
# Start defaults to node n, which is assumed to be the sink node
function relabeling_bfs(C::SparseMatrixCSC,F::SparseMatrixCSC,flowtol::Union{Float64,Int}=1e-6,start::Int64=0)

    if flowtol >= .1
        println("flowtol is a tolerance parameter for rounding small residual capacity edges to zero, and should be much smaller than $flowtol. Changing it to default value 1e-6")
        flowtol = 1e-6
    end

    # To avoid subtraction cancellation errors that may have ocurred when pushing
    # flow, when computing a bfs, round edges to zero if they are under
    # a certain tolerance
    Cf = C-F
    Cf = Cf.*(Cf.>flowtol)
    n = size(Cf,1)

    if start == 0
        start = n
    end

    rp = Cf.colptr
    ci = Cf.rowval

    N=length(rp)-1

    d = n*ones(Int64,N)
    sq=zeros(Int64,N)
    sqt=0
    sqh=0 # search queue and search queue tail/head

    # start bfs at the node "start"
    u = start
    sqt=sqt+1
    sq[sqt]=u
    d[u]=0
    while sqt-sqh>0
        sqh=sqh+1
        v=sq[sqh] # pop v off the head of the queue
        for ri=rp[v]:rp[v+1]-1
            w=ci[ri]
            if d[w] > n-1
                sqt=sqt+1
                sq[sqt]=w
                d[w]= d[v]+1
            end
        end
    end

    return d
end
