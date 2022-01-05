include("../methods.jl")
data = JSON.parsefile("fauci-email-graph.json")
##
tensor1 = JSON.parsefile("fauci-email-tensor-words.json")
R1 = temporal_reachability(data;expandstrongcomponents=false)
R = min.(R1.R,R1.R')
C1 = sort(simple_clique_heuristic(R))
@assert(data2["names"][C1] == R1.names[C1])
##
function _build_tensor_lag(tdata,subset,names;window=7,skip=6)
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
  # We want to sum over the last seven days...
  Tw = zeros(nP,nP,nT-window+1-skip,nW)
  for tw=1:nT-window+1-skip
    #@show tw+skip:(tw+skip+window-1)
    Tw[:,:,tw,:] = sum(T[:,:,tw+skip:(tw+skip+window-1),:],dims=3)
  end
  return (T=Tw, names=names[subset], dates=Date("2020-01-25") + Day.(window-1+skip:nT-1))
end
T = _build_tensor_lag(tensor1,C1,data["names"])
##
#sum(T[:,:,tw:(tw+window-1),:],dims=3)
heatmap(dropdims(sum(T.T[:,:,1,:], dims=3),dims=3))
##
using TensorDecompositions
CP(T,k) = candecomp(T,k,tuple(map(sz->randn(sz,k), size(T))...), compute_error=true, method=:ALS)
F = CP(T.T,1)
##
resids = map(k->CP(T.T,k).props[:rel_residue], 1:10)
## rank 2 seems to be a good balance... fast convergence..
F = CP(T.T,2)
##
pf = 2
plot(plot(T.dates, F.factors[3][:,pf]), plot(F.factors[1][:,pf], title="sender"),
  plot(F.factors[2][:,pf], title="receiver"), rows=2)

##
## rank 6 seems to be a good balance...
F = CP(T.T,6)
##
pf = 5
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
F = CP(T.T[2:end,2:end,:,:],10)
##
pf = 3
plot(plot(T.dates, F.factors[3][:,pf]), plot(F.factors[1][:,pf], title="sender"),
  plot(F.factors[2][:,pf], title="receiver"), rows=2)
##
@show T.names[argmax(abs.( F.factors[1][:,1]))+1] # +1 because we removed Fauci
##
@show T.names[argmax(abs.( F.factors[1][:,2]))+1] # +1 because we removed Fauci
##
for r=1:size(F.factors[1],2)
  @show r, T.names[argmax(abs.( F.factors[1][:,r]))+1] # +1 because we removed Fauci
end
##
for r=1:size(F.factors[2],2)
  @show r, T.names[argmax(abs.( F.factors[2][:,r]))+1] # +1 because we removed Fauci
end
##
Fbig = CP(T.T[2:end,2:end,:,:],20)
##
pf = 13
plot(plot(T.dates, Fbig.factors[3][:,pf]), plot(Fbig.factors[1][:,pf], title="sender"),
  plot(Fbig.factors[2][:,pf], title="receiver"), rows=2)

##
for r=1:size(Fbig.factors[1],2)
  @show r, T.names[argmax(abs.( Fbig.factors[1][:,r]))+1] # +1 because we removed Fauci
end
## Try longer windows
T = _build_tensor_lag(tensor1,C1,data["names"];window=14)
##
resids = map(k->CP(T.T[2:end,2:end,:,:],k).props[:rel_residue], 1:10)
##
F = CP(T.T[2:end,2:end,:,:],9)
##
pf = 8
plot(plot(T.dates, F.factors[3][:,pf]), scatter(F.factors[1][:,pf], title="sender"),
  scatter(F.factors[2][:,pf], title="receiver"), rows=2)
