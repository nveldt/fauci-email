## I wanted to understand how to call Gurobi directly with their API.
# This shows a little demo of calling it directly.
# It's important to use "Int32" as cind! Ahh... 

function test_gurobi(n::Int=3)

  n = 3
  obj = [-1.0,-2.0,-3.0]
  vtypes = repeat(GRB_BINARY, 3)
  aptr = Ref{Ptr{Cvoid}}()
  err = GRBnewmodel(gurobi_env, aptr, "TestFunc", 3, obj, C_NULL, C_NULL, vtypes, C_NULL)
  m = aptr[]
  try
    GRBsetintattr(m, "ModelSense", GRB_MINIMIZE)
    cind = zeros(Int32,2)
    cval = zeros(Float64,2)
    cind[1] = 0
    cind[2] = 2

    cval[1] = 1
    cval[2] = 2
    @show cind, cval
    err = GRBaddconstr(m, 2, cind, cval, GRB_EQUAL, 0.0, C_NULL)
    @show err
    GRBoptimize(m)
    GRBwrite(m, "mytest.lp");
    soln = zeros(3)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, 3, soln)
    @show soln, obj'*soln
    return soln, obj'*soln
  finally
    GRBfreemodel(m)
  end
end
test_gurobi(3)
