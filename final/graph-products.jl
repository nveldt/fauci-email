include("../methods.jl")
include("../include/Optimal_LambdaCC.jl")
include("../include/exact_conductance_jump.jl")
##
function _set_to_ind(n, S)
  ci = zeros(Int, n)
  ci[S] .= 1
  return ci
end

function _graph_products(A)
  n = size(A,1)

  c_mod = exact_modularity(A)[1]
  c_mod_zero = c_mod .- 1

  cond_S = exact_conductance(A)[1]
  c_cond = _set_to_ind(n, cond_S)

  ncut_S = exact_normalized_cut(A)[3]
  c_ncut = _set_to_ind(n, ncut_S)

  spec_S = spectral_cut(A).set
  c_spec = _set_to_ind(n, spec_S)

  return (
    modularity = c_mod_zero,
    ncut = c_ncut,
    cond = c_cond,
    spectral = c_spec,
    xy = igraph_layout(densify_graph_groups(spones!(A - Diagonal(A)), c_mod);
            random=false)
  )
end

function _compute_and_write_products(A, fn)
  p = _graph_products(A)
  open(fn, "w") do io
    JSON.print(io, p, 1)
  end
end
graphs = ["fauci-email-repliedto.json",
  "fauci-email-tofrom-5.json",
  "fauci-email-tofrom-cc-5.json",
  "fauci-email-hypergraph-projection.json",
  "fauci-email-hypergraph-projection-cc.json",]
for g in graphs
  G = _read_final(g)
  S = _simple_graph(G)
  Random.seed!(0)
  _compute_and_write_products(G.A, splitext(g)[1]*"-products-weighted.json")
  Random.seed!(0)
  _compute_and_write_products(S.A, splitext(g)[1]*"-products-simple.json")
end
  #p = _graph_products()
