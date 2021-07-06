using Gurobi
using JuMP
using SparseArrays
using LinearAlgebra
gurobi_env = Gurobi.Env()

function exact_conductance(A,outputflag = true)

    @assert(issymmetric(A))
    n = size(A,1)
    d = vec(sum(A,dims=1))
    vol = sum(d)
    I,J,V = findnz(triu(A))
    m = length(I)

    ml = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag))

    @variable(ml, x[1:n],Bin)
    @variable(ml, y[1:n])
    @variable(ml, Z[1:m],Bin)
    @variable(ml, K)

    @objective(ml, Min, K)
    @constraint(ml, sum(x) >= 1)
    @constraint(ml, sum(x[i]*d[i] for i = 1:n) - 2*sum(Z[k]*V[k] for k = 1:m) - sum(y[i]*d[i] for i = 1:n) <= 0)
    for i = 1:n
        @constraint(ml, y[i] <= x[i])
        @constraint(ml, y[i] <= K)
    end
    @constraint(ml, sum(x[i]*d[i] for i = 1:n) <= vol/2)
       
    for k = 1:m
        @constraint(ml, Z[k] <= x[I[k]])
        @constraint(ml, Z[k] <= x[J[k]])
    end

    JuMP.optimize!(ml)
    cond = JuMP.value(K)
    X = JuMP.value.(x)
    S = findall(x->x>0,X)
    Z = JuMP.value.(Z)
    Y = JuMP.value.(y)
    return S, cond, X, Z, Y
end

## test

# using MAT
# dataset = "lesmisA"
# M = matread("$dataset.mat")
# A = M["A"]

# S, cond = exact_conductance(A)
