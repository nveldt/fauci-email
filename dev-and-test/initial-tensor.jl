include("../methods.jl")
data2 = JSON.parsefile("fauci-email-graph.json")
##
tensor1 = JSON.parsefile("fauci-email-tensor-words.json")
R1 = temporal_reachability(data;expandstrongcomponents=false)
R = min.(R1.R,R1.R')
C1 = sort(simple_clique_heuristic(R))
@assert(data2["names"][C1] == R1.names[C1])
##
function _build_tensor(tdata,subset,names)
  idmap = Dict(subset .=> 1:length(subset))
  inds = tdata["indices"]
  vals = tdata["values"]
  nP = length(subset)
  nT = maximum(x->x[3], inds)+1
  nW = maximum(x->x[4], inds)+1
  T = zeros(nP,nP,nT,nW)
  for nzi in 1:length(tdata["indices"])
    (i,j,t,w) = inds[nzi].+1
    if haskey(idmap, i) && haskey(idmap, j)
      T[idmap[i],idmap[j],t,w] = vals[nzi]
    end
  end
  return (T=T, names=names[subset], dates=Date("2020-01-25") + Day.(0:nT-1))
end
T = _build_tensor(tensor1,C1,data["names"])

##
using TensorDecompositions
CP(T,k) = candecomp(T,k,tuple(map(sz->randn(sz,k), size(T))...), compute_error=true, method=:ALS)
F = CP(T.T,1)
##
resids = map(k->CP(T.T,k).props[:rel_residue], 1:10)
## rank 7 seems to be a good balance... fast convergence..
F = CP(T.T,7)
##
pf = 4
plot(plot(T.dates, F.factors[3][:,pf]), plot(F.factors[1][:,pf], title="sender"),
  plot(F.factors[2][:,pf], title="receiver"), rows=2)

##
plot(F.factors[1])
##
plot(F.factors[2])
## Fauci dominates tehse things...
# try without him...
resids = map(k->CP(T.T[2:end,2:end,:,:],k).props[:rel_residue], 1:10)

##
F = CP(T.T[2:end,2:end,:,:],8)
##
pf = 1
plot(plot(T.dates, F.factors[3][:,pf]), plot(F.factors[1][:,pf], title="sender"),
  plot(F.factors[2][:,pf], title="receiver"), rows=2)
