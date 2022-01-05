# Strongly-local code for minimizing the HCL objective
# Implemented with the thresholded linear hyperedge splitting penalty.

include("Helper_Functions.jl")
include("maxflow.jl")

"""
HyperLocal: minimizes HLC with the thresholded-linear (TL) hypergraph cut
            function, with parameter delta. In other words, the splitting function
            penalty is min { |S| , |e - S|, delta}.
            Strongly-local time

H:          Binary indicence matrix for hypergraph
Hyperedges: A list of hyperedges defining the hypergraph
order:      Order (number of nodes) in each hyperedge
d:          Degree vector, d[v] = number of hyperedges a node is in
R:          Set of nodes in seed/reference set
epsilon:    Locality parameter, must exceed vol(R)/vol(bar{R})
delta:      Threshold cut penalty.
"""
function HyperLocal(H::SparseMatrixCSC{Float64,Int64},Ht::SparseMatrixCSC{Float64,Int64},
    order::Vector{Int64},d::Vector{Float64},R::Vector{Int64},
    epsilon::Float64, delta::Float64,Rs_local::Vector{Int64},localflag::Bool=true)

    m,n = size(H)

    volA = sum(d)
    volR = sum(d[R])
    # @assert(volR <= volA/2)
    Rstrong = R[Rs_local]
    # Check Locality Parameter
    fR = volR/(volA - volR)
    @show fR, volR, volA
    if epsilon < fR
        println("Locality parameter epsilon was set too small.
        Setting it to lower bound of $fR. Computations will not be local.")
        epsilon = fR
        localflag = false
    end
    A = 0; N = 0;
    if localflag

        if volA*epsilon/volR < 10
            println("Note that vol(R)/epsilon = O(vol(G)).
            For these parameters \nit may be faster to run the algorithm
            without the locality setting.")
        end

    else
        A = tl_expansion_inc(H,order,delta)
        N = round(Int64,size(A,1))
    end

    # Store useful sets
    # Rn = hyper_neighborhood(H,R)    # get the immediate neighbors of R...
    # Rn = setdiff(Rn,R)                # ...but we exclude R itself
    Rn = get_immediate_neighbors(H,Ht,R)
    Rc = setdiff(1:n,R)               # Complement set of R
    nR = length(R)

    condR,volR, cutR = tl_cond(H,R,d,delta,volA,order)

    println("\nRunning HyperLocal")
    println("----------------------------------------")
    println("Epsilon = $epsilon \t Delta = $delta")
    println("|R| = $nR, cond(R) = $condR")
    println("-------------------------------------------------------------------------")

    S_best = R
    a_best = condR
    a_old = condR
    still_improving = true
    Iteration = 1
    while still_improving

        still_improving = false

        stepstart = time()
        if localflag
            S_new = HyperLocal_Step(H,Ht,order,R,Rn,a_best,epsilon,delta,d,Rs_local)
        else
            S_new = HLC_Step(A,R,Rc,a_best,epsilon,N,d,n,Rs_local)
        end
        stime = round(time()-stepstart,digits=1)

        a_new = hlc_tl(H,order,R,S_new,d,volA,epsilon,delta)

        if a_new < a_old
            still_improving = true
            S_best = S_new
            nS = length(S_best)
            a_old = a_new
            a_best = a_new
            println("Iter $Iteration: |S| = $nS, lcond(S) = $a_new, min-cut took $stime seconds")
        else
            println("Iter $Iteration: Algorithm converged. Last min-cut took $stime sec")
            println("-------------------------------------------------------------------------")
        end
        Iteration += 1
    end

    return S_best, a_best
end


# A non-local version of the min-cut code that works by calling the same
# subroutine, but on the entire graph all at once
function HLC_Step(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},Rbar::Vector{Int64},
    alpha::Float64, epsilon::Float64, N::Int64, d::Vector{Float64},n::Int64,Rs_local::Vector{Int64})

        Rstrong = R[Rs_local]
        # Directly set up the flow matrix
        sVec = zeros(N)
        tVec = zeros(N)
        sVec[R] .= alpha*d[R]
        sVec[Rstrong] .= N^2
        tVec[Rbar] .= alpha*epsilon*d[Rbar]
        F = maxflow(A,sVec,tVec,0)
        Src = source_nodes_min(F)[2:end].-1
        S = intersect(1:n,Src)

        return S
end

# Strongly-local subroutine for computing a minimum s-t cut
# This uses the thresholded linear splitting function for each hyperegde
function HyperLocal_Step(H::SparseMatrixCSC{Float64,Int64},Ht::SparseMatrixCSC{Float64,Int64},
    order::Vector{Int64}, R::Vector{Int64},Rn::Vector{Int64},alpha::Float64,
    epsilon::Float64,delta::Float64,d::Vector{Float64},Rs_local::Vector{Int64})

    # Map from local node indices to global node indices
    Local2Global = [R; Rn]

    n = length(d)

    # Keep track of which nodes are in the local hypergraph L
    inL = zeros(Bool,n)
    inL[Local2Global] .= true

    # Number of nodes in the local graph
    Lsize = length(Local2Global)

    # Complete nodes = nodes whose hyperedge set in the local hypergraph
    #                   is the same as their global hyperedge set
    # Incomplete nodes = everything else in the local hypergraph
    #                   (must be a neighbor of a complete node)
    #
    # Initialize the complete set to be R
    # Incomplete set is R-complement
    C_global = R
    I_global = Rn

    # Indices, in the local graph, of complete and incomplete nodes
    C_local = collect(1:length(R))
    I_local = collect(length(R)+1:Lsize)
    R_local = collect(1:length(R))
    Rstrong_local = R_local[Rs_local]

    # Get the set of hyperedges to expand around.
    # At first this is every hyperedge that touches
    # a node from R.
    Hc = H[:,C_global]
    rp_c = Hc.rowval
    # ci_c = Hc.colptr
    L_edges = unique(rp_c)

    # Binary indicence matrix for the local hypergraph (without terminal edges)
    HL = H[L_edges,Local2Global]
    order_L = order[L_edges]

    # Expand into a directed graph
    A_L = tl_expansion_inc(HL,order_L,delta)
    N_L = size(A_L,1)           # includes auxiliary nodes
    n_L = length(Local2Global)  # number of non-auxiliary nodes in A_L

    # Find the first mincut, which can be done by calling HLC_Step
    # with localized objects
    S_local = HLC_Step(A_L,C_local,I_local,alpha,epsilon,N_L,d[Local2Global],n_L,Rstrong_local)

    # Find nodes to "expand" around:
    #   any nodes in the cut set tha are "incomplete" still
    E_local = intersect(S_local,I_local)
    E_global = Local2Global[E_local]

    # ne = length(E_global)
    # println("There are $ne new nodes to expand on")

    # As long as we have new nodes to expand around, we haven't yet found
    # the global minimum s-t cut, so we continue.
    while length(E_local) > 0

        # Update which nodes are complete and which are incomplete
        C_local = [C_local; E_local]
        C_global = Local2Global[C_local]

        # Take these away from I_local
        I_local = setdiff(I_local,E_local)

        # This is better
        Nbs_of_E = get_immediate_neighbors(H,Ht,E_global)
        Lnew = setdiff(Nbs_of_E,Local2Global)
        numNew = length(Lnew)
        # Update the set of indices in L
        Local2Global = [Local2Global; Lnew]

        # Store local indices for new nodes added to L
        Lnew_local = collect((Lsize+1):(Lsize+numNew))
        Lsize = length(Local2Global)

        # These are going to be "incomplete" nodes
        I_local = [I_local; Lnew_local]
        I_global = Local2Global[I_local]

        # Now we have a new set of complete and incomplete edges,
        # we do the same thing over again to find a localize min-cut
        Hc = H[:,C_global]
        rp_c = Hc.rowval
        # ci_c = Hc.colptr
        L_edges = unique(rp_c)

        # Binary indicence matrix for the local hypergraph (without terminal edges)
        HL = H[L_edges,Local2Global]
        order_L = order[L_edges]

        # Expand into a directed graph
        A_L = tl_expansion_inc(HL,order_L,delta)
        N_L = size(A_L,1)           # includes auxiliary nodes
        n_L = length(Local2Global)  # number of non-auxiliary nodes in A_L

        # Find the first mincut, which can be done by calling HLC_Step
        # with localized objects
        R_bar_l = setdiff(1:n_L,R_local)
        S_local = HLC_Step(A_L,R_local,R_bar_l,alpha,epsilon,N_L,d[Local2Global],n_L,Rstrong_local)

        # Find nodes to "expand" around:
        #   any nodes in the cut set tha are "incomplete" still
        E_local = intersect(S_local,I_local)
        E_global = Local2Global[E_local]
        # ne = length(E_global)
        # println("There are $ne new nodes to expand on")
    end

    return Local2Global[S_local]
end
