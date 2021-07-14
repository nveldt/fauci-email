include("../methods.jl")

## Get graph
kf = false
ms = 5
mindeg = 0
G = _build_email_tofrom_graph(data; keepcc = true,
    maxset=ms, keepfauci=kf,mindegree = mindeg, emailweight = false) |> igraph_layout
drawgraph(G,markersize = 5)

##
include("../include/Optimal_LambdaCC.jl")

# Make this unweighted, if desired
A = Float64.(G.A)
fill!(A.nzval, 1)

# Exact normalized cut and modularity
# May take a while...
@time rval1= LazyExactLambdaCCGurobi(A, 0.01, true, true)

##
@time rval2= LazyExactLambdaCC(A, 0.01, true, true)
