## Wanted to look into the disconnected graph...
function _build_email_hypergraph_projection_disconnected(data;
    maxset::Int=typemax(Int), mindegree::Int=0, keepfauci=true,
    hyperedgeparts=("sender","recipients","cc"),
    emailweight::Bool=false)
  edges = Tuple{Int,Int}[]
  weights = Float64[]
  emails=data["emails"]
  names=data["names"]
  idfauci=findfirst(names .== "fauci, anthony")-1 # offset
  for thread in emails
    for email in thread
      people = Set{Int}()
      for group in hyperedgeparts
        map(p->push!(people, p), email[group])
      end
      # project the hyperedge
      if length(people) <= maxset
        for (pi,pj) in Iterators.product(people, people)
          if pi > pj
            if keepfauci || (keepfauci == false && pj!=idfauci && pi!=idfauci)
              push!(edges, (pi+1,pj+1))
              push!(weights, 1/binomial(length(people),2))
            end
          end
        end
      end
    end
  end
  A = sparse(first.(edges), last.(edges), emailweight ? weights : 1, length(names), length(names))
  A = dropzeros!(max.(A,A'))
  println("Raw data has ")
  println("  $(size(A,1)) nodes and $(nnz(A)÷2) edges")
  return (A=A, names=names)
end
G = _build_email_hypergraph_projection_disconnected(data;hyperedgeparts=("sender","recipients"))
Acc,filt = largest_component(G.A)
##
filtnodes = findall(filt.==0)
filtnames = G.names[filtnodes]
##
println.(filtnames);
## find the thread with kevin fennelly...
fid = nodeid(G, "fenne") # assumes G is set above.
for tid=1:length(data["emails"])
  for eid=1:length(data["emails"][tid])
    if data["emails"][tid][eid]["sender"] ==fid
      print_thread(data, tid)
      @show tid
      break
    end
  end
end
