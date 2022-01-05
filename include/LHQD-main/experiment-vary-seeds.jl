using Distributed
addprocs(10);

@everywhere include("evaluation.jl")
using MAT
using Random
using Profile
using JLD,FileIO

include("common.jl")
include("local-hyper.jl")

ratios = [0.1,0.08,0.04,0.02,0.01,0.008,0.004,0.002,0.001]
H,clusters = read_dataset("amazon")
records = Dict()
for ratio in ratios
    records[string(ratio)] = Dict()
    for label in keys(clusters)
        records[string(ratio)][string(label)] = Dict()
    end
end

ntrials = 10
G = LH.graph(H,10.0)
edges_to_keep = findall(G.order.<=20)
Hr = H[edges_to_keep,:]
Gr = LH.graph(Hr,1.0)
for (i,ratio) in enumerate(ratios)
    kappa = 0.025*ratios[i]
    x_eps=1.0e-16
    aux_eps=1.0e-16
    for label in keys(clusters)
        T = clusters[label][2]
        records[string(ratio)][string(label)]["LH-2.0"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=x_eps,aux_eps=aux_eps,kappa=kappa,q=2.0,ratio=ratio)
        records[string(ratio)][string(label)]["LH-1.4"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=x_eps,aux_eps=aux_eps,kappa=kappa,q=1.4,ratio=ratio)
        records[string(ratio)][string(label)]["OneHop+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_flow,epsilon=0.05,ratio=ratio,delta=G.delta,expand_num=1000)
        records[string(ratio)][string(label)]["flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_flow,epsilon=0.1,ratio=ratio,delta=G.delta,expand_ratio=0.0)
        records[string(ratio)][string(label)]["star+acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl,kappa=kappa,ratio=ratio)
        records[string(ratio)][string(label)]["clique+acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl_implicit,kappa=kappa,ratio=ratio)
        records[string(ratio)][string(label)]["weighted_clique+acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl_implicit,kappa=kappa,ratio=ratio,weighted=true)
        records[string(ratio)][string(label)]["HGCRD"] = bulk_eval_parallel(Gr,T,label,ntrials,eval_hgcrd,ratio=ratio,iterations=6)
        JLD.save("results/amazon_records_vary_seeds.jld",records)
    end
end