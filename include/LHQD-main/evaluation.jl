using Distributed
include("common.jl")
include("local-hyper.jl")
include("PageRank.jl")
include("hyperlocal_code/HyperLocal.jl")
include("SLQ.jl")
include("implicit-acl-clique-fast.jl")
include("hgcrd.jl")
include("qdsfm-ppr.jl")

"""
T:        true cluster
seed:     seed used to select random seed set
G:        HyperGraphAndDegrees
"""
function eval_lh(G,T,seed;ratio=0.01,delta=0.0,max_iters=1000000,
        x_eps=1.0e-8,aux_eps=1.0e-8,kappa=0.0025,gamma=0.1,rho=0.5,q=2.0)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    L = LH.loss_type(q,delta)
    s1 = time()
    x,r,iter = LH.lh_diffusion(G,seedset,gamma,kappa,rho,L,max_iters=max_iters,x_eps=x_eps,aux_eps=aux_eps)
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_flow(G,T,seed;ratio=0.01,epsilon=0.1,delta=1.0,expand_num=2000)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    # grownum = round(Int64,min(nT*expand_ratio,size(G.H,2)))
    OneHop = get_immediate_neighbors(G.H,G.Ht,seedset)
    grownum = round(Int64,min(expand_num,length(OneHop)))
    Rmore = BestNeighbors(G.H,G.deg,seedset,OneHop,grownum)
    R = union(Rmore,seedset)
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    s1 = time()
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,epsilon,delta,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_acl(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,
        Ga=nothing,deg=nothing,acl_delta=nothing)
    if Ga === nothing
        Ga,deg = hypergraph_to_bipartite(G)
    end
    if acl_delta === nothing
        acl_delta = G.delta
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    x = PageRank.acl_diffusion(Ga,deg,seedset,gamma,kappa)
    x ./= deg
    x = x[1:size(G.H,2)]
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,acl_delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end


function eval_acl_implicit(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,acl_delta=nothing,weighted=false,use_original=false)
    if acl_delta === nothing
        acl_delta = G.delta
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    if weighted
        Ge = ImplicitClique.ImplicitWeightedCliqueExpansion(G.H)
    else
        Ge = ImplicitClique.ImplicitCliqueExpansion(G.H)
    end
    alpha = 1/(1+gamma)
    if use_original
        x,xr = ImplicitClique.orig_acliter(Ge, alpha, seedset, kappa)
    else
        x,xr = ImplicitClique.acliter(Ge, alpha, seedset, kappa)
    end
    x ./= Ge.d
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,acl_delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end


function eval_acl_clique(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,
    Ga=nothing,deg=nothing,acl_delta=nothing,weighted=false,binary=false)
    if Ga === nothing
        Ga,deg = CliqueExpansion(G,weighted=weighted,binary=binary)
    end
    if acl_delta === nothing
        acl_delta = G.delta
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    x = PageRank.acl_diffusion(Ga,deg,seedset,gamma,kappa)
    x ./= deg
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,acl_delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_slq(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,q=1.4,delta=0.0,rho=0.5,max_iters=1000000,
        Ga=nothing,deg=nothing)
    if Ga === nothing
        Ga,deg = hypergraph_to_bipartite(G)
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    L = SLQ.loss_type(q,delta)
    x,_,_ = SLQ.slq_diffusion(SLQ.GraphAndDegrees(Ga,deg),seedset,gamma,kappa,rho,L,max_iters=max_iters)
    x = x[1:size(G.H,2)]
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_slq_clique(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,q=1.4,delta=0.0,rho=0.5,max_iters=1000000,
    Ga=nothing,deg=nothing)
    if Ga === nothing
        Ga,deg = CliqueExpansion(G)
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    L = SLQ.loss_type(q,delta)
    x,_,_ = SLQ.slq_diffusion(SLQ.GraphAndDegrees(Ga,deg),seedset,gamma,kappa,rho,L,max_iters=max_iters)
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_slq_and_flow(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,q=1.4,delta=0.0,rho=0.5,max_iters=1000000,
    Ga=nothing,deg=nothing,flow_eps=0.1)
    if Ga === nothing
        Ga,deg = hypergraph_to_bipartite(G)
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    L = SLQ.loss_type(q,delta)
    x,_,_ = SLQ.slq_diffusion(SLQ.GraphAndDegrees(Ga,deg),seedset,gamma,kappa,rho,L,max_iters=max_iters)
    x = x[1:size(G.H,2)]
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    R = cluster
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,flow_eps,1.0,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_acl_and_flow(G,T,seed;ratio=0.01,kappa=0.00025,gamma=0.1,flow_eps=0.1,
        Ga=nothing,deg=nothing,acl_delta=nothing)
    if Ga === nothing
        Ga,deg = hypergraph_to_bipartite(G)
    end
    if acl_delta === nothing
        acl_delta = G.delta
    end
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    x = PageRank.acl_diffusion(Ga,deg,seedset,gamma,kappa)
    x ./= deg
    x = x[1:size(G.H,2)]
    cond,cluster = hyper_sweepcut(G.H,x,G.deg,acl_delta,0.0,G.order,nseeds=seednum)
    R = cluster
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,flow_eps,1.0,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_lh_and_flow(G,T,seed;ratio=0.01,flow_delta=1.0,max_iters=1000000,
    x_eps=1.0e-8,aux_eps=1.0e-8,kappa=0.00025,gamma=0.1,rho=0.5,flow_eps=0.1,q=2.0)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    L = LH.loss_type(q,0.0)
    s1 = time()
    x,r,iter = LH.lh_diffusion(G,seedset,gamma,kappa,rho,L,max_iters=max_iters,x_eps=x_eps,aux_eps=aux_eps)
    cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    R = cluster
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,flow_eps,flow_delta,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_BN(G,T,seed;ratio=0.01)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),1)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    kS = nT-length(seedset)
    OneHop = get_immediate_neighbors(G.H,G.Ht,seedset)
    volA = sum(G.deg)
    s1 = time()
    B = BestNeighbors(G.H,G.deg,seedset,OneHop,kS)
    pr, re, f1 = PRF(T,B)
    cond, vol, cut = tl_cond(G.H,B,G.deg,1.0,volA,G.order)
    s2 = time()-s1
    return B,cond,cond,pr,re,f1,s2
end

function eval_BN_and_flow(G,T,seed;ratio=0.01,flow_eps=0.1)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    kS = nT-length(seedset)
    OneHop = get_immediate_neighbors(G.H,G.Ht,seedset)
    volA = sum(G.deg)
    s1 = time()
    B = BestNeighbors(G.H,G.deg,seedset,OneHop,kS)
    R = B
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,flow_eps,1.0,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,volA,G.order)
    s2 = time()-s1
    return cluster,cond,cond,pr,re,f1,s2
end

function eval_TN(G,T,seed;ratio=0.01)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    kS = nT-length(seedset)
    OneHop = get_immediate_neighbors(G.H,G.Ht,seedset)
    volA = sum(G.deg)
    s1 = time()
    B = TopNeighbors(G.H,seedset,OneHop,kS)
    pr, re, f1 = PRF(T,B)
    cond, vol, cut = tl_cond(G.H,B,G.deg,1.0,volA,G.order)
    s2 = time()-s1
    return B,cond,cond,pr,re,f1,s2
end

function eval_TN_and_flow(G,T,seed;ratio=0.01,flow_eps=0.1)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    kS = nT-length(seedset)
    OneHop = get_immediate_neighbors(G.H,G.Ht,seedset)
    volA = sum(G.deg)
    s1 = time()
    B = TopNeighbors(G.H,seedset,OneHop,kS)
    R = B
    Rs = findall(x->in(x,seedset),R)     # Force seed nodes to be in output set
    cluster, lcond = HyperLocal(G.H,G.Ht,G.order,G.deg,R,flow_eps,1.0,Rs,true)
    pr, re, f1 = PRF(T,cluster)
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,volA,G.order)
    s2 = time()-s1
    return cluster,cond,cond,pr,re,f1,s2
end

function eval_hgcrd(G,T,seed;ratio=0.01,U=3,h=3,w=2,iterations=10,alpha=1,tau=0.5)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    volA = sum(G.deg)
    s1 = time()
    lcond, cluster = capacityReleasingDiffusion(G.H,G.order,seedset,U,h,w,iterations,alpha,tau,volA,G.deg)
    pr, re, f1 = 0,0,0
    try
        pr, re, f1 = PRF(T,cluster)
    catch
        cluster = seedset
        pr, re, f1 = PRF(T,cluster)
    end
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,volA,G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end

function eval_global(G,T,seed;ratio=0.01,beta=0.05,niter=2000,h=0.01)
    nT = length(T)
    seednum = max(round(Int64,ratio*nT),5)
    Random.seed!(seed)
    p = randperm(nT)
    seedset = T[p[1:seednum]]
    s1 = time()
    x = QDSFMPageRank.qdsfmpr_ppr_euler(G.Ht, seedset, beta, niter, h)
    lcond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order,nseeds=seednum)
    pr, re, f1 = PRF(T,cluster)
    s2 = time()-s1
    cond, vol, cut = tl_cond(G.H,cluster,G.deg,1.0,sum(G.deg),G.order)
    return cluster,cond,lcond,pr,re,f1,s2
end


function make_jobs(ntrials,jobs)
    for i = 1:ntrials
        put!(jobs,i)
    end
    for i in 1:length(workers())
        put!(jobs,-1)
    end
end

function worker(G,T,jobs,results,func;kwargs...)
    while true
        seed = take!(jobs)
        if seed == -1
            break
        end
        cluster,cond,lcond,pr,re,f1,t = func(G,T,seed;kwargs...)
        put!(results,(seed,cond,lcond,pr,re,f1,t))
    end
end

function bulk_eval_parallel(G,T,label,ntrials,func;kwargs...)
    jobs = RemoteChannel(()->Channel{Int64}(ntrials+length(workers())))
    results = RemoteChannel(()->Channel{Tuple{Int64,Float64,Float64,Float64,Float64,Float64,Float64}}(ntrials))
    records = []
    make_jobs(ntrials,jobs)
    for p in workers()
        remote_do(worker,p,G,T,jobs,results,func;kwargs...)
    end
    while ntrials > 0
        ret = take!(results)
        push!(records,ret)
        ntrials -= 1
        @show label,func,ret
    end
    return records
end