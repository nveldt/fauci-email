include("methods.jl")
tensor1 = JSON.parsefile("fauci-email-tensor-words.json")
R1 = temporal_reachability(data;expandstrongcomponents=false);
R = min.(R1.R,R1.R')

# indices for people in the temportal strongly connected component
C1 = sort(simple_clique_heuristic(R))
@assert(data["names"][C1] == R1.names[C1])

## Build a tensor representation 
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

## Few emails are sent in the first week. We will restrict the time window.
datedeg = vec(sum(T.T; dims=[1,2,4]))
@show datedeg[1:6]
keepdates = vec(7:102)

## Find all the words that actually appear among these 77 people after the first week
worddeg = vec(sum(T.T[:,:,7:end,:]; dims=[1,2,3]))
keepwords = findall(x->x>0,worddeg)

# Construct a new tensor by iterating through original tensor, and only keeping entries 
# that correpond to people from C1, dates in keepdates, and words in keepwords
inds = tensor1["indices"]
vals = tensor1["values"]
names = tensor1["names"][C1]
dates = tensor1["dates"][keepdates]
words = tensor1["words"][keepwords]
idmap = Dict(C1 .=> 1:length(C1))
wordmap = Dict(keepwords .=> 1:length(keepwords))

tscindices = Vector{Any}()
tscvalues = Vector{Any}()
for nzi in 1:length(inds)
    (i,j,t,w) = inds[nzi].+1
    vl = vals[nzi]
    # check if both i & j are in temporal connected component, it's past the first week, and we have a keepword
    if haskey(idmap, i) && haskey(idmap, j) && t > 6 && haskey(wordmap,w)
        thekey = [idmap[i]-1;idmap[j]-1;t-7;wordmap[w]-1]
        push!(tscindices,thekey)
        push!(tscvalues,vl)
    end
end

## Save
description = """
 Tensor index (i, j, k, l) with value v corresponds to
 i sent 
 j an email 
 k days after 2020-01-31 and the
 lth word was closest to the embedding of the email body and
 v is the inverse of the number of recipients
 
 The third index covers 96 consecutive days
 
 There can be duplicate indices corresponding to emails on the same day
 
 CCed people are not included

 Only people in the temporal strong component are included
"""
open("fauci-email-tensor-words-tsc.json","w") do f
    JSON.print(f, Dict("indices"=>tscindices, "values"=>tscvalues, "names"=>names,"dates"=>dates,"words"=>words, "description"=>description))
end

## Load it
tsc = JSON.parsefile("fauci-email-tensor-words-tsc.json")
