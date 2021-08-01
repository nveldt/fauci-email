##
include("../methods.jl")
include("../methods_hypergraph.jl")

H = _read_final_hypergraph("fauci-email-hypergraph.json")
##

s = nodeid(H,"collins");
t = nodeid(H,"conrad");
order = vec(sum(H.H,dims = 2))
m,n = size(H.H)

## Compute a sequence of s-t cuts at the boundary between s and t nodes
Deltas = collect(1.0:0.1:5.0)
C = zeros(n,length(Deltas))     # source-sets
for i = 1:length(Deltas)
    delta = Deltas[i]

    # There are a few different choices of weight here...
    S_hyp = hyperstcut(H.H, s,t,H.weights./(order .- 1),delta);
    # S_hyp = hyperstcut(H.H, s,t,ones(m)./(order .- 1),delta);
    # S_hyp = hyperstcut(H.H, s,t,H.weights,delta);
    if in(s,S_hyp)
        S_hyp = setdiff(1:n,S_hyp)
    end
    C[S_hyp,i] .= 1
end


## Idendify nodes that change from one side to the other

Sm = vec(sum(C,dims = 2))
allx = length(Deltas)
Always0 = vec(findall(x->x==0,Sm))
Always1 = vec(findall(x->x==allx,Sm))
NeverChange = union(Always0,Always1)
Changed = setdiff(1:n,NeverChange)
I = C[Changed,:]

# Plot spy plot for just the changed nodes
a = 1; b = 2
P = a*I + b*(1 .- I)
y_label = "Boundary Nodes"
x_label = "delta"
changed_labels = H.names[Changed]
changed_labels = Vector{String}()

for k in Changed
    push!(changed_labels,split(H.names[k], ", ")[1])
end

xs = Deltas
stepx = 5
stepy = 1
ms = 4
s2 = 525
s1 = 400
spy(P,legend=:false,seriescolor=:reds,markersize=ms, size = (s2,s2),
markerstrokecolor=:white,
#bottom_margin=-20mm,
marker=:square,
xlabel=x_label, #title = "s = Collins, t = Conrad",
xtickfont=font(10), ytickfont=font(7), guidefont=font(12),titlefont=font(10),
xticks = (1:stepx:length(xs), xs[1:stepx:length(xs)]),ymirror=true,
yticks = (1:stepy:length(Changed), changed_labels[1:stepy:length(Changed)]))
##
savefig("figures/hypergraph-delta.pdf")
