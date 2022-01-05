# LHQD
code for LHQD project, see our arXiv paper 

	@Article{Liu-preprint-slq,
	  author     = {Meng Liu, Nate Veldt, Haoyu Song, Pan Li and David F. Gleich},
	  journal    = {arXiv},
	  title      = {Strongly Local Hypergraph Diffusions for Clustering and Semi-supervised Learning},
	  year       = {2020},
	  pages      = {2011.07752},
	  volume     = {cs.SI},
	  arxiv      = {http://arxiv.org/abs/2011.07752},
	  mysoftware = {https://github.com/MengLiuPurdue/LHQD},
	}


To run our code, simply `include("local-hyper.jl")` This has minimal dependencies. Then
to run the code on a small hypergraph, run

    include("common.jl")
    name = "TinyZoo"
    M = matread("small-hypergraphs/"*name*".mat")
    H = M["H"]
    good = findall(x->x > 0, order)
    H = H[good,:]
    order = vec(round.(Int64,sum(H,dims=2)))
    G = LH.graph(H,1.0) # a hypergraph object with delta=1.0 in its cut function
    q = 2.0
    L = LH.loss_type(q) # the loss-type, this is a 2-norm)
    kappa = 0.01 # value of kappa (sparsity regularization)
    gamma = 0.1 # value of gamma (regularization on seed) 
    rho = 0.5 # value of rho (KKT apprx-val)
    S = [1] $ seed set
    x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
    cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)

Datasets
-----------
The instructions on obtaining the Yelp dataset is included in the `yelp_local_algorithms` folder. We refer to reference [33] in our manuscript for instructions on obtaining the Amazon and Stack Overflow datasets. The Amazon and Stack Overflow datasets must be placed into `hypergraphs` folder.


Experiments
-----------

- Yelp experiment: check `yelp_local_algorithms` folder.
- Experiment on Amazon dataset: run `experiment-amazon.jl` to generate data, run `visualization-amazon.jl` to generate plots, run `generate-table-amazon.jl` to generate tables.
- Experiment on Stack Overflow dataset: run `experiment-stackoverflow.jl` to generate data, run `generate-table-stackoverflow.jl` to generate tables, run `visualization-stackoverflow.jl` to generate plots.
- Experiment on varying seeds: run `experiment-vary-seeds.jl` to generate data, run `visualization-vary-seeds.jl` to generate plots.
- Experiment on selecting delta: run `experiment-select-delta.jl` to generate delta, run `visualization-vary-delta.jl` to generate plots.
  



