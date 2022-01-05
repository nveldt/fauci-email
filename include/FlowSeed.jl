# This is self-contained Julia code for FlowSeed, the flow-based
# method for local clustering introduced in the paper:
#
# Flow-Based Local Graph Clustering with Better Seed Set Inclusion
# Nate Veldt, Christine Klymko, and David Gleich
# Proceedings of the 2019 SIAM International Conference on Data Mining
#
# ArXiv preprint: https://arxiv.org/abs/1811.12280
#
# The main subroutine is LocalPushRelabel.
# Unlike previous local flow methods, this repeatedly updates Phase 1 of the
# push-relabel maximum s-t flow algorithm. This phase returns a minimum s-t cut,
# which is all we need for the algorithm. The push-relabel algorithm
# is made efficient by a global relabeling heuristic.
#
# Previous flow-based methods repeatedly called a black-box min-cut
# solver and didn't use warm starts. Here we use warm starts and call a
# white-box subroutine that makes the code much faster in practice.

using SparseArrays

# This computes the precision, recall, and F1 score for a set Returned
# compared against a Target set
function PRF(Target,Returned)

    TruePos = intersect(Returned,Target)
    pr = length(TruePos)/length(Returned)
    re = length(TruePos)/length(Target)
    F1 = 2*(pr*re)/(pr+re)

    return pr, re, F1

end

# Starting from a set of seed nodes R, do a breadth first search to get
# a k-hop neighborhood of R
function neighborhood(A::SparseMatrixCSC,R::Array{Int64},k::Int64)

    rp = A.rowval
    ci = A.colptr
    n = size(A,1)

    eS = zeros(n)
    eS[R] .= 1

    # For node i, the neighbors of i are rp[ci[i]:ci[i+1]-1]
    for i = R
        neighbs = rp[ci[i]:ci[i+1]-1]
        eS[neighbs] .= 1
    end

    # This could be more efficient, but recursively calling won't take too long
    # as long as k isn't too large
    if k == 1
        return findall(x->x!=0,eS)
    else
        return neighborhood(A,findall(x->x!=0,eS),k-1)
    end

end

# For a set S in a graph with adjacency matrix A, return some information about
# S including its conductance, number of interior edges, volume, and cut.
function set_stats(A::SparseMatrixCSC{T,Int64},
    S::Vector{Int64},volA::V) where {T <: Real, V <: Real}

    if volA == 0.0
        volA = sum(A.nzval)
    end

    if length(S) == size(A,1)
        # then we have an indicator vector
        S = findall(x->x!=0,S)
        AS = A[:,S]
    else
        # then we have a subset
        @assert(minimum(S) >= 1)
        @assert(maximum(S) <= size(A,1))
        AS = A[:,S]
    end

    vol = sum(AS.nzval);
    SAS = AS[S,:]
    edges = sum(SAS.nzval);
    cut = vol-edges

    cond = cut/minimum([vol,volA-vol]);

    return cut, vol, edges, cond

end

# Compute the s-t cut score corresponding to a set S, in an augmented graph
# with source and sink node
function cutval(A::SparseMatrixCSC{Float64,Int64},S::Vector{Int64},
    R::Vector{Int64},d::Array{Float64,2},alpha::Float64,epsilon::Float64,
    volA::Float64,pR::Array{Float64},RinS::Array{Float64})

    n = size(A,1)
    if volA == 0.0
        volA = sum(A.nzval)
    end

    strongR = R[findall(x->x!=0,RinS)]
    @assert(length(setdiff(strongR,S)) == 0)    # S should contain strongR

    @assert(minimum(S) >= 1)
    @assert(maximum(S) <= size(A,1))
    AS = A[:,S];


    volS = sum(AS.nzval);
    SAS = AS[S,:]
    edges = sum(SAS.nzval);
    cutS = volS-edges

    volR = sum(d[R])

    # penalty vector, should only be nonzero for R nodes
    penalty = zeros(n)
    penalty[R] = pR.*d[R]

    RS = intersect(R,S)
    volRS = sum(d[RS])
    RnotinS = setdiff(R,RS)   # the set of nodes in R that aren't in S
    pRnotinS = sum(penalty[RnotinS])    # the penalty for excluding R nodes from A

    cutScore = cutS - alpha*volRS + alpha*volR + alpha*epsilon*(volS-volRS) + alpha*pRnotinS

    @assert(cutScore >= 0)

    relcond = cutS/(volRS - epsilon*(volS-volRS) - pRnotinS)

    return relcond
end

# The main function, which minimizes a localized variant of conductance which
# penalizes the exclusion of seed nodes from the output set.
#
# Parameters:
#
#   A = adjacency matrix for a graph
#
#   R = node indices for a seed set,
#   Rn = immediate neighbors of R
#   Rc = complement set of R
#
#   epsilon = locality parameter
#
#   pR = a length(R) vector with penalties on exluding seed nodes in R from
#       the output set. pR[i] is the penalty or excluding R[i] from the output
#
#   RinS = a length(R) zero-one vector indicating which nodes in R are stricly
#           required to be in the output set
#
#  relcondFlag = a boolean flag indicating whether to compute the relative
#                conductance score or the exact conductance score for each
#                intermediate improved set. Choosing false (i.e. updating with
#                exact conductance) will sometimes lead to fewer iterations and
#                lower conductance output, but will not actually minimize the
#                relative conductance or seed penalized conductance.
#
#  localFlag = a boolean flag indicating whether or not to use the local
#               computations. If volR is large and epsilon is small, in some
#               cases it may be better for the subroutine to perform one
#               global caluculations that multiple "local" computations.
#
# d = weighted degree vector of the graph
#
#
# volA, volR = volumes of the entire graph and seed set respectively

# FlowSeed with simplified parameters
function FlowSeed(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    relcondFlag::Bool= true,localFlag::Bool=true)

    d = sum(A,dims = 2)
    volA = sum(A.nzval)
    volR = sum(d[R])
    n = size(A,1)

    # Find one-hop neighbors of R, and get the complement set
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] .= 0
    Rc = findall(x->x!=0,inRc)             # complement of R

    if volA*epsilon/volR < 10
        localFlag = false
    end
    FlowSeed(A,R,Rn,Rc,epsilon,pR,RinS,d,volA,volR,relcondFlag,localFlag)

end

# More in depth parameters, in case one wants to run the method multiple times
# and now always recompute Rn, Rc, volA, volR, and d each time
function FlowSeed(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    Rn::Vector{Int64},Rc::Vector{Int64},epsilon::Float64,pR::Array{Float64},
    RinS::Array{Float64},d::Array{Float64},volA::Float64=0.0,volR::Float64=0.0,
    relcondFlag::Bool= true,localFlag::Bool=true)

    fR = volR/(volA - volR)
    if epsilon < fR
        println("Locality parameter epsilon was set to small. Setting it to lower bound of $fR. Computations will not be local.")
        epsilon = fR
    end

    n = size(A,1)

    if localFlag
        if volA*epsilon/volR < 10
            println("Note that vol(R)/epsilon = O(vol(G)).
            For these parameters \nit may be faster to run the algorithm
            without the locality setting.")
        end
    end

    # Call nodes that must be S the "strong seed nodes"
    localStrong = findall(x->x!=0,RinS)

    StrongSeeds = R[localStrong]
    numstrong = length(StrongSeeds)

    # If something is marked as a strong seed, put an infinite penalty
    # on excluding it from the output set
    pR[localStrong] .= Inf

    # Conductance of R
    Stats = set_stats(A,R,volA)
    alphaCurrent = Stats[4]
    # Conductance of R is same as localized seed penalized conductance of R
    # alpha2 = cutval(A,R,R,d,1.0,epsilon,volA,pR,RinS)
    # println("$alpha2, $alphaCurrent")


    println("\nEpsilon = $epsilon");
    println("There are $numstrong strong seed nodes.")
    println("The full seed set has conductance $alphaCurrent ");
    println("-------------------------------------------------------")
    BestS = R
    alph0 = 2
    alphaBest = alphaCurrent

    source = zeros(n)
    sink = zeros(n)
    dr = d[R]
    drc = d[Rc]

    while alphaCurrent < alph0

        # Prepare source-side and sink-side edge weights for the augmented
        # local flow graph
        # Seed nodes have an edge to the source of the following weight
        source[R] = alphaCurrent*(pR .+ 1).*dr

        # Non-seed nodes have an edge to the sink
        sink[Rc] = alphaCurrent*epsilon*drc

        # Compute the new min s-t cut
        if localFlag
            # Do it by repeatedly solving smaller problems, starting
            # by looking at the immediate neighbors Rn
            S = LocalPushRelabel(A,R,source,sink,Rn)
        else
            # Run a single min-cut computation on the whole graph
            S = NonLocalPushRelabel(A,R,source,sink)
        end

        if length(S) > 0 && length(S) < n

            # Check stats for new set
            if relcondFlag
                alphaS = cutval(A,S,R,d,1.0,epsilon,volA,pR,RinS)
            else
                Stats = set_stats(A,S,volA)
                alphaS = Stats[4]
            end

            if alphaS < alphaCurrent
                numS = size(S,1)
                ra = round(alphaS,digits =4)
                println("Improvement found: R-Conductance = $ra, Size = $numS")
                BestS = S
                alphaBest = alphaS
            end

        else
            alphaS = alphaCurrent
        end

        alph0 = alphaCurrent
        alphaCurrent = alphaS

    end

    SL = BestS
    sizeSL = length(SL)
    cond = alphaBest
    println("------------------------------------------------------")
    println("Final Answer: Conductance = $cond, Size = $sizeSL ")

    return SL, cond
end

# LocalPushRelabel: computes the minimumn s-t cut for a flow graph in strongly-local
#               time. It repeatedly solves localized min-cut problems.
#
# Input Parameters:
#
# A = a symmetric matrix representing an undirected graph. It can be weighted.
#
# R = a list of nodes that share an edge with the source node
#
# sWeights and tWeight store the nonnegative weight of each node to the source
# and sink. For node i, exactly one of sWeights[i] and tWeights[i] is nonzero
#
# Rn = a list of nodes not in R that neighbor a node in R
function LocalPushRelabel(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    sWeights::Array{Float64},tWeights::Array{Float64},Rn::Array{Int64})

    timer = 0.0

    n = size(A,1)
    rp = A.rowval
    ci = A.colptr

    # Now we want to locally compute maximum flows
    # C = indices of "complete" nodes in the local graph L, which are nodes
    #   whose degree in the local graph equals the degree in the global graph.
    # I = local indices of nodes that are in L, but not complete. These do
    #      share edges with one another, but only with complete nodes.

    # Initialize the complete set to be the set of nodes adjacent to the source
    C_global = R
    I_global = Rn           # everything else is incomplete
    Ac = A[C_global,:]      # set of edges from the complete set to the rest of the graph

    # We will maintain a map from indices in a local subgraph, to global indices in A.
    # These don't include the sink node in the flow graph, we are considering
    # just a growing local subgraph of A
    Local2Global = [C_global; I_global]
    # Node i in the local graph corresponds to the node with index
    # Local2Glocal[i] in the global graph A

    # Number of nodes in the local graph
    Lsize = length(Local2Global)

    # Indices, in the local graph, of complete and incomplete nodes
    C_local = collect(1:length(R))
    I_local = collect(length(R)+1:Lsize)
    numI = length(I_global)     # number of incomplete nodes

    # Build the initial local graph

    AcToI = Ac[:,I_global]     # edges between complete and incomplete nodes
    AcToc = Ac[:,C_global]     # edges between complete nodes
    L = [AcToc AcToI;
        AcToI' spzeros(numI,numI)]   # adjacency matrix for local graph

    # We distinguish between L the "local graph", and Lf, the "local flow graph"
    # which additionally contains the sink node t (as node 1).

    # In the local flow graph, each non-terminal node has either a source-side
    # or sink-side edge.
    tToL = reshape(tWeights[Local2Global],Lsize)
    sToL = reshape(sWeights[Local2Global],Lsize)

    # By adding the edges to the sink,
    # we transform the local graph L into the local flow graph Lf

    Lf = [spzeros(1,1) sparse(tToL');
         sparse(tToL) L]

    # Initialize the flow matrix; allocate space for non-zero flow values
    nLf = size(Lf,1)
    F = SparseMatrixCSC(nLf,nLf,Lf.colptr,Lf.rowval,zeros(length(Lf.rowval)))
    # Find the minimum cut for Lf.
    #
    # The first node in Lf is the sink, so offset indices of R by 1.
    start = time()
    S_local,F,excess = Main_Push_Relabel_fs(Lf,F,collect(2:length(R)+1),[0; sToL])
    timer += time()-start

    # F is a preflow that is returned. It is NOT the maximum flow for Lf.
    # S is the set of nodes in the min s-t cut of Lf. S_local are the local
    # indices in L, (not the indices in A or Lf)

    # We "expand" L around nodes in S that were previously "incomplete"
    E_local = setdiff(S_local,C_local)         # Nodes to expand around
    E_global = Local2Global[E_local]           # their global indices

    # Keep track of which nodes are in the local graph L
    inL = zeros(Bool,n)
    inL[Local2Global] .= true

    # As long as we have new nodes to expand around, we haven't yet found
    # the global minimum s-t cut, so we continue.
    while length(E_local) > 0

        # Update which nodes are complete and which are incomplete
        C_local = [C_local; E_local]
        C_global = Local2Global[C_local]

        # Take these away from I_local
        I_local = setdiff(I_local,E_local)
        I_global = Local2Global[I_local]

        # To complete nodes in E, first add all the possible edges in the
        # current local graph, so that they match the global graph edges
        # (This is one of the most expensive parts of the expansion)
        L[E_local,E_local] = A[E_global,E_global]
        L[E_local,I_local] = A[E_global,I_global]
        L[I_local,E_local] = L[E_local,I_local]'

        # Now we must expand the local graph so that NEW neighbors of E
        # are added to L
        Lnew = Vector{Int64}()
        for v = E_global
            # This extracts the neighbor list of node v from the
            # rowval and colptr vectors of the adjacency matrix
            Neighbs_of_v = rp[ci[v]:ci[v+1]-1]
            for nv = Neighbs_of_v
                if ~inL[nv]
                    inL[nv] = true
                    push!(Lnew,nv)
                end
            end
        end
        numNew = length(Lnew)

        # Store local indices for new nodes added to L
        Lnew_local = collect((Lsize+1):(Lsize+numNew))

        # These are going to be "incomplete" nodes
        I_local = [I_local; Lnew_local]

        # Expand L by adding edges from the old local graph to Lnew.
        # Note that we don't include any edges between nodes in Lnew.
        P = A[Local2Global,Lnew]
        L = [L P;
            P' spzeros(numNew,numNew)]

        # Update the set of indices in L
        Local2Global = [Local2Global; Lnew]

        # excess stores the amount of "excess" flow after a flow computation.
        #
        # Extend the excess vector to accomodate the new size of L.
        # Since Lnew were not present in the last flow computation, they
        # have zero excess.
        excess = [excess; zeros(numNew)]

        # For the next local min-cut computation, we need to know which
        # nodes come with nonzero excess. These are "active" nodes.
        ExcessNodes = findall(x->x!=0,excess)

        # Update the capacity to the sink.
        tToL = [tToL; tWeights[Lnew]]
        # Now we construct a new local flow graph, and repeat

        Lf = [spzeros(1,1) sparse(tToL');
             sparse(tToL) L]

        Fold = F    # Old flow, saved as a warm start

        # Construct an initial flow F that includes the previous flow Fold
        # as a warm start. First, we allocate space for future
        # flow.
        # (This is one of the most expensive parts of the expansion)
        nLf = size(Lf,1)

        F = SparseMatrixCSC(nLf,nLf,Lf.colptr,Lf.rowval,zeros(length(Lf.rowval)))
        F[1:Lsize+1,1:Lsize+1] = Fold

        Lsize = size(L,1)

        # Compute min s-t cut for local flow graph and see if we need to expand
        S_local,F,excess = Main_Push_Relabel_fs(Lf,F,ExcessNodes,excess)

        E_local = setdiff(S_local,C_local)     # the nodes that need completing
        E_global = Local2Global[E_local]       # their global indices

    end

    # return the global indices of the minimum cut set
    return Local2Global[S_local]

end

# A non-local version of the min-cut code that works by calling the same
# subroutine, but on the entire graph all at once
function NonLocalPushRelabel(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    sWeights::Array{Float64},tWeights::Array{Float64})

        n = size(A,1)
        # Directly set up the flow matrix
        C = [spzeros(1,1) sparse(tWeights');
            sparse(tWeights) A]

        # Allocate space for the flow we will calculate
        F = SparseMatrixCSC(n+1,n+1,C.colptr,C.rowval,zeros(length(C.rowval)))

        # R is the set of nodes with excess, and the excess
        # will come from source-side edges that are immediately saturated
        S, F, excess = Main_Push_Relabel_fs(C,F,R.+1,[0;sWeights])

        # The returned F is a preflow, not the maximum flow.
        # We are only interested in the cut.

        return S
end

# Main_Push_Relabel_fs returns a preflow F and the min s-t cut set S for the
# flow graph C. It does not solve the maximum s-t flow problem.
#
# C = the capacity matrix for the flow problem.
#   Node 1 is the sink, and there is no explicit representation of a source,
#   the preflow immediately pushes all flow from the source to create an
#   excess on nodes in the graph.
#
# F = an initial flow. It can be initialize to zero.
#
# ExcessNodes = the set of nodes which at the start begin with some positive excess
#     These can be thought of as nodes that are adjacenct to the implicit source node
#     and the edges from the source are flooded. Or they may represent nodes that
#     have a nonzero excess from the initial flow F. The indices given here
#     should account for the fact that node 1 is already reserved for the sink.
#
# excess = the vector of excess values at the start of the algorithm. If F = 0,
#   this is the vector of edge capacities from the implicit source to the graph.
#   If F != 0, then it's the excess from a previous run of the algorithm
function Main_Push_Relabel_fs(C::SparseMatrixCSC{Float64,Int64},
    F::SparseMatrixCSC{Float64,Int64},ExcessNodes::Array{Int64},excess::Array{Float64})

    # check excess node list
    # assert(countnz(excess) == length(ExcessNodes))

    # here, n includes only one terminal node, the sink
    n = size(C,1)

    height = zeros(Int64,n)      # label/height of each node
    inQ = zeros(Bool,n)          # list whether or not nodes are in the queue

    # Store adjacency list. There are ways to update this if calling
    # this function multiple times for growing local graphs, but it
    # does not appear to be a bottleneck to simply recompute frequently
    Neighbs,d = ConstructAdj(C,n)

    # We will maintain a queue of active nodes.
    Queue = Vector{Int64}()
    # An actual queue implementation is available in the DataStructures.jl
    # Julia package. The performane is nearly identical (and in some cases
    # slightly slower), thus to minimize dependency on outside packages, we
    # just use a Vector.

    # All nodes with nonzero excess are the first to be processed
    for v = ExcessNodes
        push!(Queue,v)
    end
    inQ[ExcessNodes] .= true

    # count the number of nodes that have been relabeled
    relabelings::Int64 = 0

    height = relabeling_bfs_fs(C,F)     # compute initial distance from sink

    # In the code and comments, height = distance from sink = label of node

    # Continue until the queue no longer contains any active nodes.
    while length(Queue) > 0

        u = pop!(Queue)     # Select a new active node
        inQ[u] = false      # It's no longer in the queue

        if height[u] < n    # Check that the node is still active

            # discharge flow through node u
            relabelings += discharge_fs!(C,F,Queue,u,Neighbs[u],height,excess,n,d[u],inQ)

            # if u is still active, re-place it into the queue
            if excess[u] > 0 && height[u] < n
                prepend!(Queue,u)
                inQ[u] = true
            end

        end

        # Global relabeling heuristic for push-relabel algorithm.
        # This recomputes distances between nodes and the sink
        if relabelings == n
            relabelings = 0
            dist = relabeling_bfs_fs(C,F)
            height = dist
        end

    end

    # Compute final distances from sink using BFS. Anything with distance
    # n will be the cut set.
    finalHeight = relabeling_bfs_fs(C,F)
    S = Vector{Int64}()
    for i = 2:n
        if finalHeight[i] == n
            push!(S,i-1)
        end
    end

    excess[1] = 0.0     # ignore whatever excess there was at the sink.
    return S, F, excess

end

# Discharege operation: pushes flow away from node u across admissible edges.
# If excess[u] > 0 but no admissible edges exist, we relabel u.
function discharge_fs!(C::SparseMatrixCSC{Float64,Int64},F::SparseMatrixCSC{Float64,Int64},
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    excess::Array{Float64},n::Int64,du::Int64,inQ::Array{Bool})

    vLocal::Int64 = 1
    hu = height[u]
    relabeled = 0
    while excess[u] > 0 && vLocal <= du
            v = uNeighbs[vLocal]
            if hu > height[v] && C[u,v] - F[u,v] > 0
                pushflow_fs!(C,F,Queue,u,v,excess,height,inQ,n)
                vLocal += 1
            else
                vLocal += 1
            end
    end

    if vLocal > du
        relabeled = 1
        relabel_fs!(C,F,Queue,u,uNeighbs,height,du,n)
    end

    return relabeled
end

# Relabel sets the label/height of node u to be equal to the minimum label
# such that an admissible edge exists. An edge (u,v) is admissible if
# height[u] = height[v] + 1
function relabel_fs!(C::SparseMatrixCSC{Float64,Int64},F::SparseMatrixCSC{Float64,Int64},
    Queue::Vector{Int64},u::Int64,uNeighbs::Array{Int64},height::Array{Int64},
    du::Int64,n::Int64)
   # find smallest new height making a push possible,
   # if such a push is possible at all

   min_height = Inf
   # search through the neighbors of u
   # and relabel so that height[u] = height[v] + 1 for some v in the neighborhood
   for vLocal = 1:du
       v = uNeighbs[vLocal]
       if C[u,v] - F[u,v] > 0
           min_height = min(min_height, height[v])
           height[u] = min_height + 1
       end
   end


end

# Push flow from an active node u to a node v via an admissible edge (u,v)
function pushflow_fs!(C::SparseMatrixCSC{Float64,Int64},F::SparseMatrixCSC{Float64,Int64},
    Queue::Vector{Int},u::Int64,v::Int64,excess::Array{Float64},height::Array{Int64},
    inQ::Array{Bool},n::Int64)

    send = min(excess[u], C[u,v] - F[u,v])
    F[u,v] += send
    F[v,u] -= send
    excess[u] -= send
    excess[v] += send

    # If v isn't in the queue, isn't the sink, is active, add it to the Queue
    if ~inQ[v] && v > 1 && height[v] < n
        prepend!(Queue,v)
        inQ[v] = true
    end
end

# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC{Float64,Int64},n::Int64)

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
# from each node to the sink via residual edges. Distance = n means there is no
# path to the sink. Sink node is assumed to be node 1.
function relabeling_bfs_fs(C::SparseMatrixCSC{Float64,Int64},F::SparseMatrixCSC{Float64,Int64})

    # To avoid subtraction cancellation errors that may have ocurred when pushing
    # flow, when computing a bfs we round edges to zero if they are under 1e-8
    Cf = round.((C-F),digits =6)
    n = size(Cf,1)

    rp = Cf.colptr
    ci = Cf.rowval

    N=length(rp)-1

    d = n*ones(Int64,N)
    sq=zeros(Int64,N)
    sqt=0
    sqh=0 # search queue and search queue tail/head

    # start bfs at the sink, which is node 1
    u = 1
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
