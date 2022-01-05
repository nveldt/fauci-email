using Gurobi
using JuMP
gurobi_env = Gurobi.Env()

function exact_normalized_cut(A,lam1 = 0.0,verbose = true)
    m = sum(nonzeros(A))/2
    n = size(A,1)

    if lam1 == 0.0
        lam = 1/(2*m)
        println("Starting with modularity solution.\n")
    else
        lam = lam1
        println("Starting with lambda = $lam.")
    end

    searching = true
    Clusterings = Vector{Vector{Int64}}()
    Lambdas = Vector{Float64}()
    Best_ncut = 1
    BestS = 0
    while searching

        # 1. Solve LambdaCC objective with parameter lam
        c, D =  LazyExactLambdaCC(A,lam,false,verbose)

        # Clusterings[i] = clustering obtained for parameter Lambdas[i]
        push!(Lambdas,lam)
        push!(Clusterings,c)

        # 2. Extract best scaled normalized cut, and updated lam
        lam, whichclus = compute_min_norm_cut(A,c,true)

        # 3. Save set if improvement found. Terminate otherwise
        if Best_ncut > lam
            Best_ncut = lam
            BestS = vec(findall(x->x==whichclus,c))
        else
            searching = false
        end
    end

    if length(BestS) > n/2
        BestS=setdiff(1:n, BestS)
    end

    return Clusterings, Lambdas, BestS
end

function compute_min_norm_cut(A,c,returnwhich = false)
    """
    Given a clustering c, find the minimum normalized cut score (cut(S)/[volS*vol(Sbar)]) of any cluster.
    """
    volA = sum(nonzeros(A))
    normcut = 1.0
    whichclus = 0
    for i = 1:maximum(c)
        S = vec(findall(x->x==i,c))
        cut, vol, edges, cond = set_stats(A,S,volA)
        ncut = cut/(vol*(volA-vol))
        if ncut < normcut
            normcut = ncut
            whichclus = i
        end
    end
    if returnwhich
        return normcut, whichclus
    else
        return normcut
    end
end

function exact_modularity(A)
    m = sum(nonzeros(A))/2
    d = vec(sum(A,dims = 2))
    n = size(A,1)
    lam_mod = 1/(2*m)
    c, D =  LazyExactLambdaCC(A,lam_mod,false)
    modularity_score = compute_modularity(A,d,c)
    return c, modularity_score
end


function compute_modularity(A,d,c)
    obj = 0
    m = sum(nonzeros(A)/2)
    n = size(A,1)
    for i = 1:n
        for j = 1:n
            if c[i] == c[j]
                obj += A[i,j] - d[i]*d[j]/(2*m)
            end
        end
    end

    return 1/(2*m)*obj
end

# Find exact solution to the LambdaCC objective, the degree weighted version
function LamCC_dw(A,lam)

    n = size(A,1)
    d = sum(A,dims=1)
    # Create a model that will use Gurobi to solve
    # m = Model(solver=GurobiSolver(OutputFlag = 0))

    m = Model(with_optimizer(Gurobi.Optimizer))
    # @variable(m, x[1:n,1:n], Bin)
    @variable(m, x[1:n,1:n],binary = true)

    @objective(m, Min, sum((A[i,j]-d[i]*d[j]*lam)*x[i,j] for i=1:n-1 for j = i+1:n))

    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                    @constraint(m,x[i,j] - x[i,k] - x[j,k] <= 0)
                    @constraint(m,-x[i,j] + x[i,k] - x[j,k] <= 0)
                    @constraint(m,-x[i,j] - x[i,k] + x[j,k] <= 0)
            end
        end
    end

    optimize!(m)


    # Return clustering and objective value
    D = JuMP.value.(x)
    return extractClustering(D), D
end


# extractClustering
# Given a 0,1 indicator matrix for node clustering, extract the n x 1
# cluster indicator vector
function extractClustering(x)
    # Assuming the triangle inequality results work, we just have to go through
    # each row of x, and if an element hasn't been placed in a cluster yet, it
    # it starts its own new one and adds everything that is with it.
    n = size(x,1)
    NotClustered = fill(true,n)
    c = zeros(n)
    clusnum = 1
    for i = 1:n
        if NotClustered[i]
            for j = i:n
                if x[i,j] < .01 # must be zero, they are clustered together
                    c[j] = clusnum
                    NotClustered[j] = false;
                end
            end
            clusnum +=1
        end
    end
    return round.(Int64,c)
end


# find_violations
# Given a candidate distance matrix D, iterate through all 3-tuples of nodes
# and return the tuples where triangle inequality constraints have been violated.
#
# Output is stored in vector 'violations'
#
# Heuristic speed ups:
#       - Only grab Dij, Dik, Djk once per tuple
#       - if a - b - c > 0, then a>b and a>c. By checking these first
#           we rule out a lot of situations where constraints do not need to be
#           be added, so in all test cases this gave a speed up
#
# Note that we want D to be lower triangular here, though if it is symmetric
# this will work fine as well. We just need to make sure the distance information
# in D is not just stored in the upper triangular portion
function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}})
  n = size(D,1)

  # We only need this satisfied to within a given tolerance, since the
  # optimization software will only solve it to within a certain tolerance
  # anyways. This can be tweaked if necessary.
  epsi = 1e-8
  @inbounds for i = 1:n-2
       for j = i+1:n-1
          a = D[j,i]
           for k = j+1:n
              b = D[k,i]
              c = D[k,j]
        if a - b > epsi && a - c > epsi && a-b-c > epsi
            push!(violations, (i,j,k))
                # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
        end
        if b - a > epsi && b - c > epsi && b-a-c > epsi
            push!(violations, (i,k,j))
            # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
        end

        if c - a > epsi && c-b>epsi && c-a-b > epsi
            push!(violations, (j,k,i))
            # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
        end
      end
    end
  end
end


function LazyExactLambdaCC(A,lam,outputflag = true,verbose = true)

    n = size(A,1)
    d = vec(sum(A,dims=1))
    W = zeros(n,n)
    for i = 1:n
        for j = 1:n
            W[i,j] = A[i,j]-d[i]*d[j]*lam
        end
    end
    # Create a model that will use Gurobi to solve
    # m = Model(solver=GurobiSolver(OutputFlag = 1,TimeLimit = time_limit, FeasibilityTol = FeasTol, Method = SolverMethod, Crossover = CrossoverStrategy, LogFile = OutputFile))
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag))

    @variable(m, x[1:n,1:n],Bin)

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((W[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))

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

    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
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
    for i = 1:n-1
        for j = i+1:n
            if W[i,j] > 0
                LPbound += W[i,j]*D[i,j]
            else
                LPbound += -W[i,j]*(1-D[i,j])
            end
        end
    end
    return extractClustering(D), D
end

# Use Gurobi to solve the Lambda CC objective exactly
# See Gurobi Parameter explanations online at: http://www.gurobi.com/documentation/8.0/refman/parameters.html
function LazyExactLambdaCC(A,lam,outputflag = true,verbose = true)

    n = size(A,1)
    d = vec(sum(A,dims=1))
    W = zeros(n,n)
    for i = 1:n
        for j = 1:n
            W[i,j] = A[i,j]-d[i]*d[j]*lam
        end
    end
    # Create a model that will use Gurobi to solve
    # m = Model(solver=GurobiSolver(OutputFlag = 1,TimeLimit = time_limit, FeasibilityTol = FeasTol, Method = SolverMethod, Crossover = CrossoverStrategy, LogFile = OutputFile))
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag))

    @variable(m, x[1:n,1:n],Bin)

    # Minimize the weight of disagreements (relaxation)
    @objective(m, Min, sum((W[i,j])*x[i,j] for i=1:n-1 for j = i+1:n))

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

    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    @constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                end
            end
        end
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
    for i = 1:n-1
        for j = i+1:n
            if W[i,j] > 0
                LPbound += W[i,j]*D[i,j]
            else
                LPbound += -W[i,j]*(1-D[i,j])
            end
        end
    end
    return extractClustering(D), D
end


function LazyExactLambdaCCGurobi(A,lam,outputflag = true,verbose = true)
    n = size(A,1)
    d = vec(sum(A,dims=1))
    W = zeros(n,n)
    for i = 1:n
        for j = 1:n
            W[i,j] = A[i,j]-d[i]*d[j]*lam
        end
    end
    #nvars = div(n*(n-1), 2)
    # build varmap and obj
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    for j=1:n
        for i=1:j-1
            nvars += 1
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    vtypes = repeat(GRB_BINARY, nvars)
    aptr = Ref{Ptr{Cvoid}}()
    err = GRBnewmodel(gurobi_env, aptr, "LazyLambdaCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    m = aptr[]

    try

    cind = Int32[0,0,0]
    cval = Float64[0,0,0]
    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    #@constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                    cind[1] = vmap[min(j,k),max(j,k)]
                    cind[2] = vmap[min(i,k),max(i,k)]
                    cind[3] = vmap[min(i,j),max(i,j)]
                    cval[1] = 1
                    cval[2] = -1
                    cval[3] = -1
                    error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                end
            end
        end
    end

    # Find intial first solution
    if verbose
        println("First round of optimization")
    end
    #JuMP.optimize!(m)
    GRBoptimize(m)
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
    obj = robj[]
    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
    D = zeros(n,n)
    for j=1:n
        for i=1:j-1
            D[i,j] = soln[vmap[i,j]+1]
            D[j,i] = soln[vmap[i,j]+1]
        end
    end

    while true
        # x will naturally be upper triangular, but 'find_violations'  wants lower triangular
          #D = Matrix(JuMP.value.(x)')

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
             cind[1] = vmap[v[1],v[2]]
             cind[2] = vmap[min(v[1],v[3]),max(v[1],v[3])]
             cind[3] = vmap[min(v[2],v[3]),max(v[2],v[3])]
             cval[1] = 1
             cval[2] = -1
             cval[3] = -1
             #@show cind, cval
             error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
             #z = sparsevec(cind.+1, cval, nvars)
             #@show z'*soln
             #@constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
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
         GRBupdatemodel(m)
         err = GRBoptimize(m)
         #@show err
         robj = Ref{Float64}(0.0)
         GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
         obj = robj[]
         #@show "obj after solve", obj
         GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
         for j=1:n
             for i=1:j-1
                 D[i,j] = soln[vmap[i,j]+1]
                 D[j,i] = soln[vmap[i,j]+1]
             end
         end

     end
     return extractClustering(D), D
 finally
     #@show "freemodel here!"
     GRBfreemodel(m)
 end
end
