## Nate asked why we have self-loops,
# so I wanted to look at why...
include("../methods.jl")
##
for thread in data["emails"]
  edges,weights = _build_tofrom_edges(thread)
  if any(first.(edges) .== last.(edges))
    print_thread(thread, data["names"])
  end
end
