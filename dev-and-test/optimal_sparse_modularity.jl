
##
include("../include/Optimal_LambdaCC.jl")
##
function badtris(A::SparseMatrixCSC)
    nz = A.nzval
    cp = A.colptr
    rv = A.rowval
    n = size(A,1)
    tris = Tuple{Int,Int,Int}[]
    neighind = zeros(Int,n)
    for i=1:n
        for nzi=cp[i]:cp[i+1]-1
            j = rv[nzi]
            if j > i
                neighind[j] = nzi
            end
        end

        for nzi=cp[i]:cp[i+1]-1
            j = rv[nzi]
            if j > i
                for nzi2 = cp[j]:cp[j+1]-1
                    k = rv[nzi2]
                    if k > j && neighind[k] > 0 # (i->j->k->i) by symmetry
                        parity = nz[nzi]*nz[nzi2]*nz[neighind[k]]
                        if parity < 0
                            push!(tris, (i,j,k))
                        end
                    end
                end
            end
        end

        # reset indicator
        for nzi=cp[i]:cp[i+1]-1
            j = rv[nzi]
            neighind[j] = 0
        end
    end
    println("$(length(tris)) bad triangles")
    return tris
end



function LazyExactSparseModularity(B::SparseMatrixCSC,outputflag = true,verbose = true)

    W = copy(B)
    n = size(W,1)

    edgevals = collect(zip(findnz(M)[1:3]...))

    # Create a model that will use Gurobi to solve
    # m = Model(solver=GurobiSolver(OutputFlag = 1,TimeLimit = time_limit, FeasibilityTol = FeasTol, Method = SolverMethod, Crossover = CrossoverStrategy, LogFile = OutputFile))
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag))

    @variable(m, x[1:n,1:n],Bin)

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((v)*x[i,j] for (i,j,v) in edgevals))

    # Constraints: 0 <= x_ij <= 1 for all node pairs (i,j)
    # for i = 1:n-1
    #     for j = i:n
    #         @constraint(m,x[i,j] <= 1)
    #         @constraint(m,x[i,j] >= 0)
    #     end
    # end

    # If we solve with no triangle inequality constraints, we would get the
    # solution D[i,j] = D[j,i] = 1 if A[i,j] = 0 and D[i,j] = D[j,i] = 1
    # otherwise.
    #
    # If we then checked all 3-tuples of nodes, we'd find that the only time
    # triangle constraints were violated is at bad triangles (++- triangles,
    # a tuple that is an unclosed wedge).
    # Hence, from the very beginning we know we must include
    # a triangle inequality constraint at every bad triangle.
    println()
    println("Adding bad triangles...")
    btris = badtris(W)
    for (i,j,k) in btris
        @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
    end

    # Find intial first solution
    if verbose
        println("First round of optimization")
    end
    JuMP.optimize!(m)

    while true
        # x will naturally be upper triangular, but 'find_violations'  wants lower triangular
          D = Matrix(JuMP.value.(x)')

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()

         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)
         if verbose
            print("Adding in $numvi violated constraints.")
         end

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             @constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end
         if numvi == 0
             if verbose
                println(" Optimal solution found.")
             end
             break
         end
         if verbose
            println(" And re-solving the ILP.")
         end
         JuMP.optimize!(m)
     end

    # Return the objective score and the distance matrix
    D = JuMP.value.(x)
    LPbound = 0.0
    for (i,j,v) in zip(findnz(W)[1:3]...)
        LPbound += v*D[i,j]
    end
    return extractClustering(D), D, LPbound
end


function exact_temporal_modularity(B)
    #m = sum(nonzeros(A))/2
    c, D,bound =  LazyExactSparseModularity(B)
end

Tcc_ids = vec([ 672 728 729 730 733 734 826 735 737 751 792 818 133 80 81 83 88 102 123 124 125 79 137 140 159 161 162 220 229 34 2 3 4 5 7 22 30 32 239 35 36 48 54 67 68 77 536 551 556 559 576 584 508 613 635 639 666 1 677 694 350 241 245 276 278 292 293 294 299 368 375 399 412 480 498 503 795])
T = build_temporal_graphs(data; subset=Tcc_ids)
M, slices = expand_to_slicetime_matrix(T)
c,score = exact_temporal_modularity(M)
