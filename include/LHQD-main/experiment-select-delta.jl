using Distributed
using JLD,FileIO
addprocs(10);

@everywhere include("evaluation.jl")

include("common.jl")
include("local-hyper.jl")

H,clusters = read_dataset("amazon",min_size=10000,max_size=40000)
ntrials = 30
deltas = [1.0,5.0,10.0,50.0,100.0,500.0,1000.0,5000.0,10000.0]

labels = [29,4,21,10,22]

records = Dict{String,Dict{String,Dict{String,Vector{Tuple}}}}()
for label in labels
    records[string(label)] = Dict{String,Dict{String,Vector{Tuple}}}()
    for delta in deltas
        records[string(label)][string(delta)] = Dict{String,Vector{Tuple}}()
    end
end

for label in labels
    for delta in deltas
        T = clusters[label][2]
        G = LH.graph(H,delta)
        kappa = 0.00025
        ratio = 0.01
        records[string(label)][string(delta)]["LH-2.0"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=1.0e-10,aux_eps=1.0e-16,kappa=kappa,q=2.0,ratio=ratio)
        @everywhere GC.gc()
    end
    JLD.save("results/amazon_records_vary_delta.jld",records)
end


H,clusters = read_dataset("stackoverflow")
ntrials = 30
deltas = [1.0,5.0,10.0,50.0,100.0,500.0,1000.0,5000.0,10000.0]

records = Dict{String,Dict{String,Dict{String,Vector{Tuple}}}}()
for label in collect(keys(records))[1:5]
    records[string(label)] = Dict{String,Dict{String,Vector{Tuple}}}()
    for delta in deltas
        records[string(label)][string(delta)] = Dict{String,Vector{Tuple}}()
    end
end

for label in collect(keys(records))[1:5]
    for delta in deltas
        T = clusters[label][2]
        G = LH.graph(H,delta)
        kappa = 0.0025
        ratio = 0.02
        records[string(label)][string(delta)]["LH-2.0"] = bulk_eval_parallel(G,T,label,ntrials,eval_lh,x_eps=1.0e-10,aux_eps=1.0e-16,kappa=kappa,q=2.0)
        @everywhere GC.gc()
    end
    JLD.save("results/stackoverflow_records_vary_delta.jld",records)
end