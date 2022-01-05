using MutableArithmetics
using Random
using MAT
using SparseArrays
using DataStructures

function PRF(Target,Returned)

    if length(Returned) == 0
        pr = 0; re = 0; F1 = 0
    else
        TruePos = intersect(Returned,Target)
        pr = length(TruePos)/length(Returned)
        re = length(TruePos)/length(Target)
        F1 = 2*(pr*re)/(pr+re)

        if length(TruePos) == 0
            F1 = 0
        end
    end

    return pr, re, F1

end

function capacityReleasingDiffusion(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, refNodes::Vector{Int64}, U::Int64, h::Int64, w::Int64, iterations::Int64, alpha::Int64, tau::Float64, volA::Float64, d::Vector{Float64})
	n = size(d)[1]
	sumDegree = sum(d)
	delta = spzeros(Float64, n, 1)
	refDegree = 0
	dh = Dict{Int64, Float64}()
	for refNode in refNodes
		if !haskey(dh, refNode)
			e = findnz(H[:,refNode])[1]
			dh[refNode] = sum(order[e])
		end
		delta[refNode] = w*dh[refNode]
		refDegree += dh[refNode]
	end
	
	condBest = 100
	cutBest = Vector{Float64}
	for i in 1:iterations
		println("iteration = ", i)
		fv = spzeros(Float64, n, 1)
		ex = spzeros(Float64, n, 1)
		l = spzeros(Int64, n, 1)
		unitFlowHyperedge(H, order, delta, U, h, w, alpha, fv, ex, l, dh)

		cutTemp, condTemp = roundUnitFlow(H, order, l, volA, d)
		if condTemp < condBest
			condBest = condTemp
			cutBest = cutTemp
		end
		excess = 0
		fvsnz = findnz(fv)[1]
		for fvnz in fvsnz
			excess += ex[fvnz]
			delta[fvnz] = w*(fv[fvnz]-ex[fvnz])
		end

		stopCondition = refDegree*(w ^i)*tau
		if excess > stopCondition 
        		break
		end
        	sum_ = 0
        	for fvnz in fvsnz
			e = findnz(H[:,fvnz])[1]
        		dv = sum(order[e])
			if fv[fvnz] >= dv
				sum_ += dv 
			end
       	end
		#if sum_ > sumDegree/3 
        	#	break
        	#end
    	end
	return condBest, cutBest
end



function tl_cond(H::SparseMatrixCSC,S::Vector{Int64},order::Vector{Int64})

    volS = 0
    # Check the cut
    HS = H[:,S]
    sumHS = sum(HS,dims = 2)  # Count number of S nodes in each hyperedge
    inds = findall(x->x>0,sumHS)    # Extract hyperedges with > 0 nodes from S
    ES = sumHS[inds]
    verts = order[inds]               # Get the size of these hyperedges
    cutval = 0.0
    for j = 1:length(ES)
        if ES[j] != verts[j]
                cutval += 1
        end
	   volS += ES[j]
    end

    return cutval/min(volS, sum(order)-volS), volS, cutval
end

function roundUnitFlow(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, l::SparseMatrixCSC{Int64,Int64}, volA::Float64, d::Vector{Float64})
	n = size(H)[2]
	lsnz = findnz(l)[1]
	labels = spzeros(Int64, n, n)
	keys = Vector{Int64}()
	for lnz in lsnz 
		labels[l[lnz], labels[l[lnz],1]+2] = lnz
		labels[l[lnz],1] = labels[l[lnz],1] + 1
		append!(keys, l[lnz])
	end
    
	keys = sort(unique(keys), rev=true)
	S = Vector{Int64}()
	cond = 100
	SFinal = Vector{Int64}()
	for level in keys
    		if labels[level,1] == 0
    			continue
    		else
			for i=2:labels[level,1]+1
				append!(S, labels[level, i])
			end
			temp_cond, vol, cut = tl_cond(H,S,order)
			if temp_cond < cond
				cond = temp_cond
				SFinal =  deepcopy(S)
			end 
    		end
    	end
	return SFinal, cond
end


function pushHyperedge(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, f::Dict{String,Float64}, fv::SparseMatrixCSC{Float64,Int64}, v::Int64, u::Array{Int64,1}, l::SparseMatrixCSC{Int64,Int64}, ex::SparseMatrixCSC{Float64,Int64}, U::Int64, n::Int64, w::Int64, alpha::Int64, idx::String, motifSize::Int64,dh::Dict{Int64, Float64})
	if ex[v] < (motifSize - 1)
		return false
	end
	pushed = false
	if !haskey(f, idx)
		f[idx] = 0
	end
	r = U - f[idx]
	countDownhill = 0
	for i in 1:(motifSize - 1)	
		if l[v] > l[(u[i])]
			countDownhill = countDownhill+1
			if countDownhill >= alpha
				break
			end
		end
	end
 	if (r > 0) && countDownhill >= alpha
				psi = min(floor(ex[v]/(motifSize-1)),r)
				for i in 1:(motifSize - 1)
					if !haskey(dh, u[i])
						e = findnz(H[:,u[i]])[1]
						dh[u[i]] = sum(order[e])
					end
					psi = min(psi, w*dh[u[i]]-fv[(u[i])])
				end
				if psi > 0
					f[idx] += psi;
					for i = 1:(motifSize - 1)
						fv[v] -= psi;
						fv[(u[i])] += psi
					end
					pushed = true
				end
	end
	return pushed
end

function pushRelabelHyperedge(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, f::Dict{String,Float64}, fv::SparseMatrixCSC{Float64,Int64}, U::Int64, v::Int64, currentv::SparseMatrixCSC{Int64,Int64}, ex::SparseMatrixCSC{Float64,Int64}, l::SparseMatrixCSC{Int64,Int64}, n::Int64, w::Int64, alpha::Int64, dh::Dict{Int64, Float64})
	index = currentv[v]
	hv = findnz(H[:,v])[1]
	if index + 1 > length(hv)
		l[v] = l[v] + 1
		relabelled = true
		currentv[v] = 0
		return false, true, []
	end
	e = H[hv[index+1],:]
	u = findnz(e)[1]
	u = setdiff(u, v)
	pushed = pushHyperedge(H,order,f,fv,v,u,l,ex,U,n,w,alpha,string(v, ",", index),length(u)+1, dh)
	relabelled = false
	if (!pushed) 
		if index+1 < length(hv)
			currentv[v] = currentv[v] + 1
		else 
			l[v] = l[v] + 1
			relabelled = true
			currentv[v] = 0
		end
	end
	return pushed,relabelled,u
end

function updateExcessHyperEdge(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, fv::SparseMatrixCSC{Float64,Int64}, v::Int64, ex::SparseMatrixCSC{Float64,Int64}, motifSize::Int64, dh::Dict{Int64, Float64})
	if !haskey(dh, v)
		e = findnz(H[:,v])[1]
		dh[v] = sum(order[e])
	end
	ex[v] = max(fv[v]-dh[v],0.0)
end


function unitFlowHyperedge(H::SparseMatrixCSC{Float64,Int64}, order::Vector{Int64}, Delta::SparseMatrixCSC{Float64,Int64}, U::Int64, h::Int64, w::Int64, alpha::Int64, fv::SparseMatrixCSC{Float64,Int64}, ex::SparseMatrixCSC{Float64,Int64}, l::SparseMatrixCSC{Int64,Int64}, dh::Dict{Int64, Float64})
	Q = PriorityQueue() 
	n = size(H)[2]
	currentv =  spzeros(Int64, n, 1)
	f =  Dict{String,Float64}()
	deltanz = findnz(Delta)[1]
	for delta in deltanz
		if !haskey(dh, delta)
			e = findnz(H[:,delta])[1]
			dh[delta] = sum(order[e])
		end
		fv[delta] = Delta[delta]
		degreeVal = dh[delta]
		if Delta[delta] > degreeVal 
			l[delta] = 1
			enqueue!(Q, delta, l[delta])
			ex[delta] = Delta[delta] - degreeVal
		end
	end
	cnt = 0
	while !isempty(Q)
		cnt = cnt + 1
		if cnt % 10000 == 0
			println("cnt = ", cnt, " , ", length(Q))
		end
		v = peek(Q)[1]
		tmpResult = pushRelabelHyperedge(H,order,f,fv,U,v,currentv,ex,l,n,w,alpha,dh)
		u = tmpResult[3]
		if tmpResult[1]
				updateExcessHyperEdge(H,order,fv,v,ex, length(u)+1, dh)
				if ex[v] == 0
					dequeue!(Q,v)
				end
				for i=1:length(u)
					updateExcessHyperEdge(H,order,fv,u[i],ex, length(u)+1, dh)
					if ex[u[i]] > 0 && !haskey(Q, u[i]) && l[u[i]] < h
						enqueue!(Q, u[i], l[u[i]])
					end
				end
		end
		if tmpResult[2] 
			if l[v] < h
				dequeue!(Q,v)
				enqueue!(Q, v, l[v])
			else
				dequeue!(Q,v)
			end
		end
	end
end





