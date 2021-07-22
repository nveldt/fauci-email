using Distributed
addprocs(3);

@everywhere include("evaluation.jl")
using MAT
using Random
using Profile
using JLD,FileIO
using Statistics

include("common.jl")
include("local-hyper.jl")

using PyCall
@pyimport collections as py_collections
H,clusters = read_dataset("amazon")
G = LH.graph(H,1.0)
kappa = 0.00025
ratio = 0.01
flow_eps = 0.05
ntrials = 30
x_eps=1.0e-10
aux_eps=1.0e-16
expand_num=1000
edges_to_keep = findall(G.order.<=20)
Hr = H[edges_to_keep,:]
Gr = LH.graph(Hr,1.0)

records = Dict()
for label in labels
    records[string(label)] = Dict()
end

for label in keys(clusters)
    T = clusters[label][2]
    records[string(label)]["clique+acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl_implicit,kappa=kappa,ratio=ratio)
    records[string(label)]["weighted_clique+acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl_implicit,kappa=kappa,ratio=ratio,weighted=true)
    records[string(label)]["acl"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl,kappa=kappa,ratio=ratio)
    records[string(label)]["LH-2.0"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=x_eps,aux_eps=aux_eps,kappa=kappa,q=2.0,ratio=ratio)
    records[string(label)]["LH-1.4"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=x_eps,aux_eps=aux_eps,kappa=kappa,q=1.4,ratio=ratio)
    records[string(label)]["LH-2.0+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh_and_flow,x_eps=x_eps,
        aux_eps=aux_eps,kappa=kappa,q=2.0,flow_eps=flow_eps,ratio=ratio)
    records[string(label)]["HGCRD"] = bulk_eval_parallel(Gr,T,label,ntrials,eval_hgcrd,ratio=ratio,iterations=6)
    records[string(label)]["BN"] = bulk_eval_parallel(G,T,label,ntrials,eval_BN,ratio=ratio)
    records[string(label)]["TN"] = bulk_eval_parallel(G,T,label,ntrials,eval_TN,ratio=ratio)
    records[string(label)]["OneHop+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_flow,ratio=ratio,delta=G.delta,expand_num=expand_num,epsilon=flow_eps)
    acl_flow_med_old = median([t[end-1] for t in records[string(label)]["acl+flow"]])
    records[string(label)]["acl+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_acl_and_flow,
        kappa=kappa,ratio=ratio,flow_eps=flow_eps)
    acl_flow_med_new = median([t[end-1] for t in records[string(label)]["acl+flow"]])
    records[string(label)]["BN+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_BN_and_flow,
        ratio=ratio,flow_eps=flow_eps)
    records[string(label)]["TN+flow"] = bulk_eval_parallel(G,T,label,ntrials,eval_TN_and_flow,
        ratio=ratio,flow_eps=flow_eps)
    @everywhere GC.gc()
    save("results/amazon_records.jld",records)
end