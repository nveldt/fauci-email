using LinearAlgebra
using SparseArrays

module LH

using SparseArrays, DataStructures, LinearAlgebra, ProgressMeter

"""
H:          |E| x |V| binary incidence matrix for the hypergraph
Ht:         A list of hyperedges defining the hypergraph, |V| x |E|
order:      Order (number of nodes) in each hyperedge
d:          Degree vector, d[v] = number of hyperedges a node is in
delta:      Threshold cut penalty.
"""
struct HyperGraphAndDegrees{
        T<: Union{Float32,Float64,Int32,Int64},
        Ti <: Union{Int,Int32,Int64}}   # T is the type of edges,
  H::SparseMatrixCSC{T,Ti}
  Ht::SparseMatrixCSC{T,Ti}
  delta::T
  deg::Vector{T}
  order::Vector{Ti}
end

abstract type EdgeLoss{T} end

struct QHuberLoss{T} <: EdgeLoss{T}
    q::T
    delta::T
end

struct TwoNormLoss{T} <: EdgeLoss{T}
end

""" This function isn't type stable, so don't use it except it outer codes. """
function loss_type(q::T, delta::T) where T
    if q == 2.0
        return TwoNormLoss{T}()
    else
        return QHuberLoss{T}(q, delta)
    end
end

loss_type(q) = loss_type(q,0.0)

minval(f, L::QHuberLoss) = f^(1/(L.q-1))
minval(f, L::TwoNormLoss) = sqrt(f)

function loss_gradient(x::T, L::QHuberLoss{T}) where T
    if abs(x) < L.delta
        return L.delta^(L.q-2)*max(x,0)
    else
        return max(x,0)^(L.q-1)
    end
end

function loss_gradient(x::T, L::TwoNormLoss{T}) where T
    return max(x,0)
end

function loss_gradient_rev(y::T, L::QHuberLoss{T}) where T
    ret = y^(1/(L.q-1))
    if ret < L.delta
        return y/(L.delta^(L.q-2))
    else
        return ret
    end
end

function loss_gradient_rev(y::T, L::TwoNormLoss{T}) where T
    return y
end

function loss_function(x::T, L::QHuberLoss{T}) where T
    if abs(x) < L.delta
        return 0.5*(L.delta^(L.q-2))*(x^2)
    else
        return (abs(x)^L.q)/L.q+(0.5-1/L.q)*(L.delta^L.q)
    end
end

function loss_function(x::T, L::TwoNormLoss{T}) where T
    return 0.5*(x^2)
end

function graph(H::SparseMatrixCSC,delta)
    Ht = sparse(H')
    order = vec(round.(Int64,sum(H,dims=2)))
    d = vec(sum(H,dims=1))
    return HyperGraphAndDegrees(H,Ht,delta,d,order)
end

"""
buf_a, buf_b, buf_a_val, buf_b_val will be ordered such that all the a or b 
that are smaller than x is on the left
i.e. if x = 3, buf_a = [4,2,1,5,6], then buf_a is reordered as [2,1,6,5,4]
"""
function _buffer_neighbors!(xi::T, a::Vector, b::Vector, H::SparseMatrixCSC,
        i::Int, buf_a::Vector{T}, buf_b::Vector{T}, buf_a_vals::Vector{T}, buf_b_vals::Vector{T}) where T
    nneighs = H.colptr[i+1]-H.colptr[i]
    min_val,max_val = 1.0,0.0
    a_thd,b_thd = 1.0,1.0
    sum_a,sum_b,sum_weight = 0.0,0.0,0.0
    a_left_ptr,b_left_ptr = 1,1
    a_right_ptr,b_right_ptr = nneighs,nneighs
    for (iter,k) in enumerate(H.colptr[i]:(H.colptr[i+1]-1))
        j = H.rowval[k]
        val = H.nzval[k]
        min_val = min(min_val,a[j])
        max_val = max(max_val,a[j])
        if a[j] <= xi
            buf_a[a_left_ptr] = a[j]
            buf_a_vals[a_left_ptr] = T(val)
            sum_a += a[j]*T(val)
            sum_weight += T(val)
            a_left_ptr += 1
        else
            buf_a[a_right_ptr] = a[j]
            buf_a_vals[a_right_ptr] = T(val)
            a_thd = min(a_thd,a[j])
            a_right_ptr -= 1
        end
        min_val = min(min_val,b[j])
        max_val = max(max_val,b[j])
        if b[j] <= xi
            buf_b[b_left_ptr] = b[j]
            buf_b_vals[b_left_ptr] = T(val)
            b_left_ptr += 1
        else
            buf_b[b_right_ptr] = b[j]
            buf_b_vals[b_right_ptr] = T(val)
            sum_b += b[j]*T(val)
            sum_weight += T(val)
            b_thd = min(b_thd,b[j])
            b_right_ptr -= 1
        end
    end
    # @assert(a_left_ptr-a_right_ptr==1)
    # @assert(b_left_ptr-b_right_ptr==1)
    return nneighs,a_left_ptr-1,b_left_ptr-1,min_val,max_val,a_thd,b_thd,sum_a,sum_b,sum_weight
end

"""
values in buf_x is ordered such that from left to right, x that is smaller than both ai and bi
appears first, then appears x that is smaller than ai but larger than bi and then x that is larger
than both ai and bi
"""
function _buffer_neighbors!(x::Vector, Ht::SparseMatrixCSC,
        i::Int, buf_x::Vector{T}, buf_x_vals::Vector{T}, ai::T, bi::T) where T
    # @assert(ai >= bi)
    nneighs = Ht.colptr[i+1]-Ht.colptr[i]
    min_val,max_val = 1.0,0.0
    left_ptr,right_ptr = 1,nneighs
    sum_x_a,sum_x_b = 0.0,0.0
    sum_weight_a,sum_weight_b = 0.0,0.0
    x_a_thd,x_b_thd = 1.0,1.0
    tmp = []
    for (iter,k) in enumerate(Ht.colptr[i]:(Ht.colptr[i+1]-1))
        j = Ht.rowval[k]
        val = Ht.nzval[k]
        if x[j] <= bi
            buf_x[left_ptr] = x[j]
            buf_x_vals[iter] = T(val)
            sum_x_b += x[j]*T(val)
            sum_weight_b += T(val)
            left_ptr += 1
        elseif x[j] > ai
            buf_x[right_ptr] = x[j]
            buf_x_vals[iter] = T(val)
            sum_x_a += x[j]*T(val)
            sum_weight_a += T(val)
            x_b_thd = min(x[j],x_b_thd)
            x_a_thd = min(x[j],x_a_thd)
            right_ptr -= 1
        else
            x_b_thd = min(x[j],x_b_thd)
            push!(tmp,k)
        end
        min_val = min(min_val,x[j])
        max_val = max(max_val,x[j])
    end
    b_left_ptr = left_ptr
    for k in tmp
        j = Ht.rowval[k]
        val = Ht.nzval[k]
        buf_x[left_ptr] = x[j]
        buf_x_vals[left_ptr] = T(val)
        left_ptr += 1
    end
    a_left_ptr = left_ptr
    # @assert(left_ptr-right_ptr==1)
    return nneighs,a_left_ptr-1,b_left_ptr-1,min_val,max_val,sum_x_a,sum_x_b,sum_weight_a,sum_weight_b,x_a_thd,x_b_thd
end

"""
reorder values in val_buf using thd as piviot
"""
function _reorder_buffer(val_buf::AbstractVector{T},weight_buf::AbstractVector{T},
        thd,si,ei,dir,curr_ptr) where T
    sum_val_changed = 0
    sum_weight_changed = 0
    if dir == 1
        for i = si:ei
            if val_buf[i] < thd
                sum_val_changed -= val_buf[i]*weight_buf[i]
                sum_weight_changed -= weight_buf[i]
                val_buf[curr_ptr],val_buf[i] = val_buf[i],val_buf[curr_ptr]
                weight_buf[curr_ptr],weight_buf[i] = weight_buf[i],weight_buf[curr_ptr]
                curr_ptr += 1
            end
        end
        curr_ptr -= 1
    else
        for i = si:-1:ei
            if val_buf[i] > thd
                sum_val_changed += val_buf[i]*weight_buf[i]
                sum_weight_changed += weight_buf[i]
                val_buf[curr_ptr],val_buf[i] = val_buf[i],val_buf[curr_ptr]
                weight_buf[curr_ptr],weight_buf[i] = weight_buf[i],weight_buf[curr_ptr]
                curr_ptr -= 1
            end
        end
    end
    return curr_ptr,sum_val_changed,sum_weight_changed
end

function _eval_residual_i!(xi::T, di::T, dx::T, seed::Bool,
        neigh_a::AbstractVector{T}, neigh_b::AbstractVector{T},
        neigh_a_vals::AbstractVector{T}, neigh_b_vals::AbstractVector{T},
        L::QHuberLoss{T}, gamma::T, left_ptr, right_ptr,direction) where T
    num_a_small,num_b_small = left_ptr
    num_a_boundary,num_b_boundary = right_ptr
    if direction > 0
        curr_ptr = num_a_small+1
        si,ei = curr_ptr,num_a_boundary
    else
        curr_ptr = num_a_small
        si,ei = curr_ptr,num_a_boundary+1
    end
    num_a_small,_,_ = _reorder_buffer(
        neigh_a,neigh_a_vals,xi+dx,si,ei,direction,curr_ptr)
    if direction > 0
        curr_ptr = num_b_small+1
        si,ei = curr_ptr,num_b_boundary
    else
        curr_ptr = num_b_small
        si,ei = curr_ptr,num_b_boundary+1
    end
    num_b_small,_,_ = _reorder_buffer(
        neigh_b,neigh_b_vals,xi+dx,si,ei,direction,curr_ptr)
    
    ri_new = zero(T)
    for k in 1:num_a_small
      ri_new -= neigh_a_vals[k]*loss_gradient(xi+dx-neigh_a[k],L)/gamma
    end
    for k in (num_b_small+1):length(neigh_b)
      ri_new += neigh_b_vals[k]*loss_gradient(neigh_b[k]-xi-dx,L)/gamma
    end
    if seed
      ri_new += di*loss_gradient(1-xi-dx,L)
    else
      ri_new -= di*loss_gradient(xi+dx,L)
    end
    new_ptr = (num_a_small,num_b_small)
    return ri_new,new_ptr
end

# function _eval_residual_i!(xi::T, di::T, dx::T, seed::Bool,
#     neigh_a::AbstractVector{T}, neigh_b::AbstractVector{T},
#     neigh_a_vals::AbstractVector{T}, neigh_b_vals::AbstractVector{T},
#     L::QHuberLoss{T}, gamma::T) where T
#     ri_new = zero(T)
#     for k in 1:length(neigh_a)
#         ri_new -= neigh_a_vals[k]*loss_gradient(xi+dx-neigh_a[k],L)/gamma
#     end
#     for k in 1:length(neigh_b)
#         ri_new += neigh_b_vals[k]*loss_gradient(neigh_b[k]-xi-dx,L)/gamma
#     end
#     if seed
#         ri_new += di*loss_gradient(1-xi-dx,L)
#     else
#         ri_new -= di*loss_gradient(xi+dx,L)
#     end
#     return ri_new
# end

function _eval_residual_i!(xi::T, di::T, dx::T, seed::Bool,
    neigh_a::AbstractVector{T}, neigh_b::AbstractVector{T},
    neigh_a_vals::AbstractVector{T}, neigh_b_vals::AbstractVector{T},
    L::TwoNormLoss{T}, gamma::T, direction, last_ptr, boundary_ptr, sum_a, sum_b, sum_weight) where T
    nneighs = length(neigh_a)
    num_a_small,num_b_small = last_ptr
    num_a_boundary,num_b_boundary = boundary_ptr
    if direction > 0
        curr_ptr = num_a_small+1
        si,ei = curr_ptr,num_a_boundary
    else
        curr_ptr = num_a_small
        si,ei = curr_ptr,num_a_boundary+1
    end
    num_a_small,sum_val_changed,sum_weight_changed = _reorder_buffer(
        neigh_a,neigh_a_vals,xi+dx,si,ei,direction,curr_ptr)
    sum_a -= sum_val_changed
    sum_weight -= sum_weight_changed
    if direction > 0
        curr_ptr = num_b_small+1
        si,ei = curr_ptr,num_b_boundary
    else
        curr_ptr = num_b_small
        si,ei = curr_ptr,num_b_boundary+1
    end
    num_b_small,sum_val_changed,sum_weight_changed = _reorder_buffer(
        neigh_b,neigh_b_vals,xi+dx,si,ei,direction,curr_ptr)
    sum_b += sum_val_changed
    sum_weight += sum_weight_changed
    ri_new = sum_a/gamma+sum_b/gamma-(sum_weight/gamma)*(xi+dx)-di*(xi+dx)+seed*di
    new_ptr = (num_a_small,num_b_small)
    # needs to re-consider a_large and b_large
    # if direction == 1
    #     a_left_ptr = num_a_small+1
    #     b_left_ptr = num_b_small+1
    #     for i = (num_a_small+1):num_a_boundary
    #         if neigh_a[i] < xi+dx
    #             sum_a += neigh_a_vals[i]*neigh_a[i]
    #             sum_weight += neigh_a_vals[i]
    #             neigh_a[a_left_ptr],neigh_a[i] = neigh_a[i],neigh_a[a_left_ptr]
    #             neigh_a_vals[a_left_ptr],neigh_a_vals[i] = neigh_a_vals[i],neigh_a_vals[a_left_ptr]
    #             a_left_ptr += 1
    #         end
    #     end
    #     for i = (num_b_small+1):num_b_boundary
    #         if neigh_b[i] < xi+dx
    #             sum_b -= neigh_b_vals[i]*neigh_b[i]
    #             sum_weight -= neigh_b_vals[i]
    #             neigh_b[b_left_ptr],neigh_b[i] = neigh_b[i],neigh_b[b_left_ptr]
    #             neigh_b_vals[b_left_ptr],neigh_b_vals[i] = neigh_b_vals[i],neigh_b_vals[b_left_ptr]
    #             b_left_ptr += 1
    #         end
    #     end
    #     new_ptr = (a_left_ptr-1,b_left_ptr-1)
    # # needs to re-consider a_small and b_small
    # else
    #     a_right_ptr = num_a_small
    #     b_right_ptr = num_b_small
    #     for i = num_a_small:-1:(num_a_boundary+1)
    #         if neigh_a[i] > xi+dx
    #             sum_a -= neigh_a_vals[i]*neigh_a[i]
    #             sum_weight -= neigh_a_vals[i]
    #             neigh_a[a_right_ptr],neigh_a[i] = neigh_a[i],neigh_a[a_right_ptr]
    #             neigh_a_vals[a_right_ptr],neigh_a_vals[i] = neigh_a_vals[i],neigh_a_vals[a_right_ptr]
    #             a_right_ptr -= 1
    #         end
    #     end
    #     for i = num_b_small:-1:(num_b_boundary+1)
    #         if neigh_b[i] > xi+dx
    #             sum_b += neigh_b_vals[i]*neigh_b[i]
    #             sum_weight += neigh_b_vals[i]
    #             neigh_b[b_right_ptr],neigh_b[i] = neigh_b[i],neigh_b[b_right_ptr]
    #             neigh_b_vals[b_right_ptr],neigh_b_vals[i] = neigh_b_vals[i],neigh_b_vals[b_right_ptr]
    #             b_right_ptr -= 1
    #         end
    #     end
    #     new_ptr = (a_right_ptr,b_right_ptr)
    # end
    # ri_new = sum_a/gamma+sum_b/gamma-(sum_weight/gamma)*(xi+dx)-di*(xi+dx)+seed*di
    return ri_new,new_ptr,sum_a,sum_b,sum_weight
end

"""
sum_a is the weighted sum of aj that is smaller than xi
sum_b is the weighted sum of bj that is larger than xi
sum_weight is the sum of weight for xi
a_thd is the smallest aj that is larger than xi
b_thd is the smallest bj that is larger than xi
"""
function dxi_solver!(G::HyperGraphAndDegrees,x::Vector{T},
        a::Vector{T},b::Vector{T},kappa::T,epsilon::T,gamma::T,r::Vector{T},
        seedset,rho::T,i::Int,L::TwoNormLoss{T},
        buf_a::Vector,buf_b::Vector,buf_a_vals::Vector,buf_b_vals::Vector) where T
    is_seed = i in seedset
    di = G.deg[i]
    found_dxi = false
    H = G.H
    ri_new = rho*kappa*di
    nneighs::Int,num_a_small::Int,num_b_small::Int,min_val,max_val,a_thd,b_thd,sum_a,sum_b,sum_weight = _buffer_neighbors!(
        x[i],a,b,H,i,buf_a,buf_b,buf_a_vals,buf_b_vals)
    # case 0: assume xi+dxi doesn't change the relative order of a, b, xi
    tmp1 = sum_b/gamma+sum_a/gamma-ri_new+is_seed*di
    tmp2 = sum_weight/gamma+di
    tmp3 = tmp1/tmp2
    if tmp3 <= b_thd && b_thd != 1 && tmp3 <= a_thd && a_thd != 1
        return tmp3-x[i], ri_new, a_thd, b_thd, sum_weight
    end
    # case 1: assume xi+dxi < min_val
    if num_a_small == 0 && num_b_small == 0
        tmp0 = sum(@view(buf_b[1:nneighs]).*@view(buf_b_vals[1:nneighs]))
        tmp1 = tmp0/gamma
        tmp1 -= ri_new
        tmp1 += is_seed*di
        tmp2 = di+di/gamma
        if tmp1/tmp2 <= min_val
            a_thd,b_thd = minimum(@view(buf_a[1:nneighs])),minimum(@view(buf_b[1:nneighs]))
            return tmp1/tmp2-x[i], ri_new, a_thd, b_thd, di
        end
    end
    # case 2: assume xi+dxi > max_val
    tmp0 = sum(@view(buf_a[1:nneighs]).*@view(buf_a_vals[1:nneighs]))
    tmp1 = tmp0/gamma
    tmp1 -= ri_new
    tmp1 += is_seed*di
    tmp2 = di+di/gamma
    if tmp1/tmp2 >= max_val
        return tmp1/tmp2-x[i], ri_new, 1.0, 1.0, di
    end
    # case 3: min_val < x+dx < max_val 
    left_ptr = (num_a_small,num_b_small) # at least this number of a or b will be smaller than x+dx
    right_ptr = (nneighs,nneighs) # at most this number of a or b will be smaller than x+dx
    dx_min,dx_max = max(min_val-x[i],0),max_val-x[i]
    direction = 1 # direction=1 means in the last bisect, the dx_min is evaluated, direction=-1 means dx_max is evaluated
    last_ptr = left_ptr
    boundary_ptr = right_ptr
    iters = 0
    while (found_dxi == false) && (left_ptr[1] < right_ptr[1] || left_ptr[2] < right_ptr[2])
        dx_mid = dx_max/2+dx_min/2
        ri_new,ptr_new,sum_a,sum_b,sum_weight = _eval_residual_i!(x[i], T(di), dx_mid, is_seed,
            @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
            @view(buf_b_vals[1:nneighs]),L,gamma,direction,last_ptr,boundary_ptr,sum_a,sum_b,sum_weight) 
        if ri_new < rho*kappa*di
            dx_max = dx_mid
            right_ptr = ptr_new
            boundary_ptr = left_ptr
            direction = -1
        elseif ri_new > rho*kappa*di
            dx_min = dx_mid
            left_ptr = ptr_new
            boundary_ptr = right_ptr
            direction = 1
        else
            found_dxi = true
        end
        last_ptr = ptr_new
        iters += 1
    end
    # now we know exactly which a or b will be smaller than x+dx
    # dx can then be solved immediately
    dxi = (sum_a/gamma+sum_b/gamma+is_seed*di-rho*kappa*di)/(sum_weight/gamma+di)-x[i]
    if left_ptr[1]+1 > nneighs
        a_thd = 1
    else
        a_thd = minimum(@view(buf_a[(left_ptr[1]+1):nneighs]))
    end
    if left_ptr[2]+1 > nneighs
        b_thd = 1
    else
        b_thd = minimum(@view(buf_b[(left_ptr[2]+1):nneighs]))
    end
    return dxi,rho*kappa*di,a_thd,b_thd,sum_weight
end

function dxi_solver!(G::HyperGraphAndDegrees,x::Vector{T},
        a::Vector{T},b::Vector{T},kappa::T,epsilon::T,gamma::T,r::Vector{T},
        seedset,rho::T,i::Int,L::QHuberLoss{T},
        buf_a::Vector,buf_b::Vector,buf_a_vals::Vector,buf_b_vals::Vector,thd1,thd2) where T
    # @show thd1,thd2
    di = G.deg[i]
    found_dxi = false
    H = G.H
    nneighs::Int,num_a_small::Int,num_b_small::Int,min_val,max_val = _buffer_neighbors!(
        x[i],a,b,H,i,buf_a,buf_b,buf_a_vals,buf_b_vals)

    nbisect = 0

    ri_new = r[i]
    dx_min = 0
    thd_min = min(thd1,thd2)
    thd_max = max(thd1,thd2)
    thd = thd_max
    dx = thd
    left_ptr = (num_a_small,num_b_small)
    right_ptr = (nneighs,nneighs) 
    direction = 1
    ri_new,new_ptr = _eval_residual_i!(x[i], T(di), dx, i in seedset,
        @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
        @view(buf_b_vals[1:nneighs]),L,gamma,left_ptr,right_ptr,direction)
    if ri_new < 0
        ri_new = r[i]
        thd = thd_min
        left_ptr = (num_a_small,num_b_small)
        right_ptr = (nneighs,nneighs) 
        last_ptr = left_ptr
    else
        last_ptr = left_ptr
        left_ptr = new_ptr
    end
    last_dx = 0
    # @show left_ptr,right_ptr,buf_a[1:nneighs],buf_b[1:nneighs],x[i]+dx
    ratio = 10 # 2020-05-27 switched this ratio from 2 to 10
    while ri_new > rho*kappa*di
        dx = thd
        ri_new,new_ptr = _eval_residual_i!(x[i], T(di), dx, i in seedset,
            @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
            @view(buf_b_vals[1:nneighs]),L,gamma,left_ptr,right_ptr,direction)
        # @show left_ptr,right_ptr,buf_a[1:nneighs],buf_b[1:nneighs],x[i]+dx
        last_ptr = left_ptr
        left_ptr = new_ptr
        last_dx = dx_min
        dx_min = thd
        thd *= ratio
        nbisect += 1
    end
    dx_min = last_dx
    dx_max = thd/ratio
    left_ptr = last_ptr
    right_ptr = new_ptr
    last_ptr = new_ptr
    boundary_ptr = left_ptr
    direction = -1

    dx_mid = 0
    iters = 0
    while ((found_dxi == false && dx_max - dx_min > epsilon) || (ri_new < 0)) && (dx_min<(dx_min+dx_max)/2<dx_max)
    # while ((found_dxi == false && dx_max - dx_min > epsilon) || (ri_new < 0))
        dx_mid = dx_max/2+dx_min/2
        ri_new,new_ptr = _eval_residual_i!(x[i], T(di), dx_mid, i in seedset,
            @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
            @view(buf_b_vals[1:nneighs]),L, gamma,last_ptr,boundary_ptr,direction)
        # if iters < 100
        #     @show iters,ri_new,left_ptr,rho*kappa*di,right_ptr,new_ptr,dx_min,dx_max,x[i]+dx_mid
        # end
        #@show buf_a[1:nneighs],buf_b[1:nneighs],x[i]+dx_mid,rho*kappa*di
        if ri_new < rho*kappa*di
            dx_max = dx_mid
            right_ptr = new_ptr
            boundary_ptr = left_ptr
            direction = -1
        elseif ri_new > rho*kappa*di
            dx_min = dx_mid
            left_ptr = new_ptr
            boundary_ptr = right_ptr
            direction = 1
        else
            found_dxi = true
        end
        last_ptr = new_ptr
        iters += 1
    end
    # @show "dx",iters
    if dx_mid == 0
        dxi = dx_max
    else
        dxi = dx_mid
    end
    return dxi,ri_new
end

# function dxi_solver!(G::HyperGraphAndDegrees,x::Vector{T},
#     a::Vector{T},b::Vector{T},kappa::T,epsilon::T,gamma::T,r::Vector{T},
#     seedset,rho::T,i::Int,L::QHuberLoss{T},
#     buf_a::Vector,buf_b::Vector,buf_a_vals::Vector,buf_b_vals::Vector,thd1,thd2) where T
#     di = G.deg[i]
#     found_dxi = false
#     H = G.H
#     nneighs::Int,num_a_small::Int,num_b_small::Int,min_val,max_val = _buffer_neighbors!(
#         x[i],a,b,H,i,buf_a,buf_b,buf_a_vals,buf_b_vals)

#     nbisect = 0

#     ri_new = r[i]
#     dx_min = 0
#     thd_min = min(thd1,thd2)
#     thd_max = max(thd1,thd2)
#     thd = thd_max
#     dx = thd
#     ri_new = _eval_residual_i!(x[i], T(di), dx, i in seedset,
#         @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
#         @view(buf_b_vals[1:nneighs]),L, gamma)
#     if ri_new < 0
#         ri_new = r[i]
#         thd = thd_min
#     end
#     last_dx = 0

#     ratio = 10 # 2020-05-27 switched this ratio from 2 to 10
#     while ri_new > rho*kappa*di
#         dx = thd
#         ri_new = _eval_residual_i!(x[i], T(di), dx, i in seedset,
#             @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
#             @view(buf_b_vals[1:nneighs]),L, gamma)
#         last_dx = dx_min
#         dx_min = thd
#         thd *= ratio
#         nbisect += 1
#     end
#     dx_min = last_dx
#     dx_max = thd/ratio

#     dx_mid = 0
#     iters = 0
#     init_iter = true
#     while (found_dxi == false && dx_max - dx_min > epsilon) || (ri_new < 0) || init_iter
#         init_iter = false
#         dx_mid = dx_max/2+dx_min/2
#         ri_new = _eval_residual_i!(x[i], T(di), dx_mid, i in seedset,
#             @view(buf_a[1:nneighs]),@view(buf_b[1:nneighs]), @view(buf_a_vals[1:nneighs]),
#             @view(buf_b_vals[1:nneighs]),L, gamma)
#         #@show buf_a[1:nneighs],buf_b[1:nneighs],x[i]+dx_mid,rho*kappa*di
#         if ri_new < rho*kappa*di
#             dx_max = dx_mid
#         elseif ri_new > rho*kappa*di
#             dx_min = dx_mid
#         else
#             found_dxi = true
#         end
#         iters += 1
#     end
#     return dx_mid,ri_new
# end

"""
fast_dxi_solver! assumes the relative order of x,a,b doesn't change
"""
function fast_dxi_solver!(x::Vector{T},r::Vector{T},i::Int,di::T,rho::T,kappa::T,gamma::T,
        sums_weight::Vector{T},as_thd::Vector{T},bs_thd::Vector{T}) where T
    ri_new = rho*kappa*di
    dxi = (r[i] - ri_new)/(di + sums_weight[i]/gamma)
    a_thd,b_thd = as_thd[i],bs_thd[i]
    is_invalid = a_thd == 1 || b_thd == 1 || x[i]+dxi > a_thd || x[i]+dxi > b_thd
    return dxi,ri_new,is_invalid
end

function _eval_da_db_residual!(G::HyperGraphAndDegrees,ai::T,bi::T,dai::T,neigh_x::AbstractVector{T},
    neigh_x_vals::AbstractVector{T},L::TwoNormLoss{T},last_ptr,boundary_ptr,direction,sum_x_a,sum_x_b,sum_weight_a,
    sum_weight_b) where T
    num_a_small,num_b_small = last_ptr
    num_a_boundary,num_b_boundary = boundary_ptr
    if direction > 0
        curr_ptr = num_a_small+1
        si,ei = curr_ptr,num_a_boundary
    else
        curr_ptr = num_a_small
        si,ei = curr_ptr,num_a_boundary+1
    end
    num_a_small,sum_val_changed,sum_weight_changed = _reorder_buffer(
        neigh_x,neigh_x_vals,ai+dai,si,ei,direction,curr_ptr)
    sum_x_a += sum_val_changed
    sum_weight_a += sum_weight_changed
    res_ai = sum_x_a-sum_weight_a*(ai+dai)
    dbi = ai+dai-res_ai/G.delta-bi
    if direction > 0
        curr_ptr = num_b_small+1
        si,ei = curr_ptr,num_a_small
    else
        curr_ptr = num_b_small
        si,ei = curr_ptr,num_b_boundary+1
    end
    num_b_small,sum_val_changed,sum_weight_changed = _reorder_buffer(
        neigh_x,neigh_x_vals,bi+dbi,si,ei,direction,curr_ptr)
    sum_x_b -= sum_val_changed
    sum_weight_b -= sum_weight_changed
    res_bi = -sum_x_b+sum_weight_b*(bi+dbi)
    new_ptr = (num_a_small,num_b_small)
    return new_ptr,dbi,res_ai,res_bi,sum_x_a,sum_weight_a,sum_x_b,sum_weight_b
end

function _eval_da_db_residual!(G::HyperGraphAndDegrees,ai::T,bi::T,dai::T,neigh_x::AbstractVector{T},
    neigh_x_vals::AbstractVector{T},L::QHuberLoss{T},last_ptr,boundary_ptr,direction) where T
    num_a_small,num_b_small = last_ptr
    num_a_boundary,num_b_boundary = boundary_ptr
    if direction > 0
        curr_ptr = num_a_small+1
        si,ei = curr_ptr,num_a_boundary
    else
        curr_ptr = num_a_small
        si,ei = curr_ptr,num_a_boundary+1
    end
    num_a_small,_,_ = _reorder_buffer(
        neigh_x,neigh_x_vals,ai+dai,si,ei,direction,curr_ptr)
    res_ai = 0.0
    for k in (num_a_small+1):length(neigh_x)
        res_ai += neigh_x_vals[k]*loss_gradient(neigh_x[k]-ai-dai,L)
    end
    dbi = ai+dai-loss_gradient_rev(res_ai/G.delta,L)-bi
    if direction > 0
        curr_ptr = num_b_small+1
        si,ei = curr_ptr,num_a_small
    else
        curr_ptr = num_b_small
        si,ei = curr_ptr,num_b_boundary+1
    end
    num_b_small,_,_ = _reorder_buffer(
        neigh_x,neigh_x_vals,bi+dbi,si,ei,direction,curr_ptr)
    res_bi = 0.0
    for k in 1:num_b_small
        res_bi += neigh_x_vals[k]*loss_gradient(bi+dbi-neigh_x[k],L)
    end
    new_ptr = (num_a_small,num_b_small)
    return new_ptr,dbi,res_ai,res_bi
end

# function _eval_da_db_residual!(G::HyperGraphAndDegrees,ai::T,bi::T,dai::T,neigh_x::AbstractVector{T},
#     neigh_x_vals::AbstractVector{T},L::QHuberLoss{T}) where T
#     res_ai = 0.0
#     for k in 1:length(neigh_x)
#         res_ai += neigh_x_vals[k]*loss_gradient(neigh_x[k]-ai-dai,L)
#     end
#     dbi = ai+dai-loss_gradient_rev(res_ai/G.delta,L)-bi
#     res_bi = 0.0
#     for k in 1:length(neigh_x)
#         res_bi += neigh_x_vals[k]*loss_gradient(bi+dbi-neigh_x[k],L)
#     end
#     return dbi,res_ai,res_bi
# end

function dai_dbi_solver!(G::HyperGraphAndDegrees,x::Vector{T},buf_x::Vector{T},buf_x_vals::Vector{T},
    a::Vector{T},b::Vector{T},epsilon::T,i::Int,L::TwoNormLoss{T},thd) where T
    ai,bi = a[i],b[i]
    Ht = G.Ht
    H = G.H
    nneighs::Int,num_a_small::Int,num_b_small::Int,min_val,max_val,sum_x_a,sum_x_b,sum_weight_a,sum_weight_b,x_a_thd,x_b_thd = _buffer_neighbors!(
        x,Ht,i,buf_x,buf_x_vals,ai,bi)
    # case 0: assume ai+dai and bi+dbi doesn't change the relative order of ai, bi, x
    dbi = (sum_x_a+sum_x_b+sum_weight_a*sum_x_b/G.delta)/(sum_weight_a+sum_weight_b+sum_weight_a*sum_weight_b/G.delta)-bi
    bi = bi+dbi
    ai = bi+sum_weight_b*bi/G.delta-sum_x_b/G.delta
    if ai <= x_a_thd && x_a_thd != 1 && bi <= x_b_thd && x_b_thd != 1
        return ai-a[i],dbi,x_a_thd,x_b_thd,sum_weight_a,sum_weight_b
    end
    ai,bi = a[i],b[i]
    # case 1: assume a+da < min_val => b+db < min_val, not possible
    # case 2: assume a+da > max_val, not possible either
    # case 3: min_val < a+da < max_val
    # empirically, da_max should be around c*dxi where c is a constant
    dai = thd
    da_min,da_max = 0,dai
    res_ai = 0.0
    left_ptr = (num_a_small,num_b_small)
    right_ptr = (nneighs,nneighs)
    direction = 1 # direction=1 means in the last iteration, the da_min is evaluated, direction=-1 means da_max is evaluated
    last_ptr = left_ptr
    boundary_ptr = right_ptr
    new_ptr,dbi,res_ai,res_bi,sum_x_a,sum_weight_a,sum_x_b,sum_weight_b = _eval_da_db_residual!(
        G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
        direction,sum_x_a,sum_x_b,sum_weight_a,sum_weight_b)
    ratio = 10.0
    while res_bi <= res_ai || dbi < 0
        da_min = dai
        dai *= ratio
        da_max = dai
        left_ptr = new_ptr
        last_ptr = new_ptr
        boundary_ptr = right_ptr
        direction = 1
        new_ptr,dbi,res_ai,res_bi,sum_x_a,sum_weight_a,sum_x_b,sum_weight_b = _eval_da_db_residual!(
            G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
            direction,sum_x_a,sum_x_b,sum_weight_a,sum_weight_b)
    end
    right_ptr = new_ptr
    last_ptr = new_ptr
    boundary_ptr = left_ptr
    direction = -1
    while (left_ptr[1] < right_ptr[1] || left_ptr[2] < right_ptr[2])
        dai = (da_min+da_max)/2
        new_ptr,dbi,res_ai,res_bi,sum_x_a,sum_weight_a,sum_x_b,sum_weight_b = _eval_da_db_residual!(
            G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
            direction,sum_x_a,sum_x_b,sum_weight_a,sum_weight_b)
        if res_bi > res_ai && dbi >= 0
            da_max = dai
            right_ptr = new_ptr
            boundary_ptr = left_ptr
            direction = -1
        else
            da_min = dai
            left_ptr = new_ptr
            boundary_ptr = right_ptr
            direction = 1
        end
        last_ptr = new_ptr
    end
    # now we know exactly which ai+dai or bi+dbi will be smaller than x
    # dai and dbi can then be solved immediately
    dbi = (sum_x_a+sum_x_b+sum_weight_a*sum_x_b/G.delta)/(sum_weight_a+sum_weight_b+sum_weight_a*sum_weight_b/G.delta)-bi
    bi = bi+dbi
    dai = bi+sum_weight_b*bi/G.delta-sum_x_b/G.delta-ai
    if left_ptr[1]+1 > nneighs
        x_a_thd = 1
    else
        x_a_thd = minimum(@view(buf_x[(left_ptr[1]+1):nneighs]))
    end
    if left_ptr[2]+1 > nneighs
        x_b_thd = 1
    else
        x_b_thd = minimum(@view(buf_x[(left_ptr[2]+1):nneighs]))
    end
    return dai,dbi,x_a_thd,x_b_thd,sum_weight_a,sum_weight_b
end

function dai_dbi_solver!(G::HyperGraphAndDegrees,x::Vector{T},buf_x::Vector{T},buf_x_vals::Vector{T},
        a::Vector{T},b::Vector{T},epsilon::T,i::Int,L::QHuberLoss{T},thd) where T
        ai,bi = a[i],b[i]
        Ht = G.Ht
        nneighs::Int,num_a_small::Int,num_b_small::Int,_,_,_,_,_,_,_,_ = _buffer_neighbors!(
            x,Ht,i,buf_x,buf_x_vals,ai,bi)
        dai = thd
        da_min,da_max = 0,dai
        left_ptr = (num_a_small,num_b_small)
        right_ptr = (nneighs,nneighs)
        direction = 1 # direction=1 means in the last iteration, the da_min is evaluated, direction=-1 means da_max is evaluated
        last_ptr = left_ptr
        boundary_ptr = right_ptr
        new_ptr,dbi,res_ai,res_bi = _eval_da_db_residual!(
            G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
            direction)
        ratio = 10.0
        while res_bi <= res_ai || dbi < 0
            da_min = dai
            dai *= ratio
            da_max = dai
            left_ptr = new_ptr
            last_ptr = new_ptr
            boundary_ptr = right_ptr
            direction = 1
            new_ptr,dbi,res_ai,res_bi = _eval_da_db_residual!(
                G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
                direction)        
        end
        right_ptr = new_ptr
        last_ptr = new_ptr
        boundary_ptr = left_ptr
        direction = -1
        while (da_max - da_min > epsilon) && (da_min<(da_min+da_max)/2<da_max)
            dai = (da_min+da_max)/2
            new_ptr,dbi,res_ai,res_bi = _eval_da_db_residual!(
                G,ai,bi,dai,@view(buf_x[1:nneighs]),@view(buf_x_vals[1:nneighs]),L,last_ptr,boundary_ptr,
                direction)
            if res_bi > res_ai && dbi >= 0
                da_max = dai
                right_ptr = new_ptr
                boundary_ptr = left_ptr
                direction = -1
            else
                da_min = dai
                left_ptr = new_ptr
                boundary_ptr = right_ptr
                direction = 1
            end
            last_ptr = new_ptr
        end
        # dai = thd
        # da_min,da_max = 0,dai
        # dbi,residual_ai,residual_bi = _eval_da_db_residual!(G,ai,bi,dai,@view(buf_x[1:nneighs]),
        #     @view(buf_x_vals[1:nneighs]),L)
        # ratio = 10.0
        # while residual_bi <= residual_ai || dbi < 0
        #     da_min = dai
        #     dai *= ratio
        #     da_max = dai
        #     dbi,residual_ai,residual_bi = _eval_da_db_residual!(G,ai,bi,dai,@view(buf_x[1:nneighs]),
        #     @view(buf_x_vals[1:nneighs]),L)
        # end
        # init_iter = true
        # while da_max - da_min > epsilon || init_iter
        #     init_iter = false
        #     dai = (da_min+da_max)/2
        #     dbi,residual_ai,residual_bi = _eval_da_db_residual!(G,ai,bi,dai,@view(buf_x[1:nneighs]),
        #         @view(buf_x_vals[1:nneighs]),L)
        #     if residual_bi > residual_ai && dbi >= 0
        #         da_max = dai
        #     else
        #         da_min = dai
        #     end
        # end
        return dai,max(0,dbi)
end

"""
fast_dai_dbi_solver! assumes the relative order of x,a,b doesn't change
"""
function fast_dai_dbi_solver!(i::Int,j::Int,dxi::T,x::Vector{T},a::Vector{T},b::Vector{T},sums_weight_a::Vector{T},
        sums_weight_b::Vector{T},xs_a_thd::Vector{T},xs_b_thd::Vector{T},delta::T) where T
    sum_weight_a,sum_weight_b = sums_weight_a[j],sums_weight_b[j]
    x_a_thd,x_b_thd = xs_a_thd[j],xs_b_thd[j]
    dx1,dx2 = 0,0
    if x[i]-dxi <= b[j]
        dx1 = min(dxi,b[j]-x[i]+dxi)
    end
    if x[i]-dxi > a[j]
        dx2 = min(x[i]-dxi-a[j],dxi)
    end
    dai = (dx1+dx2+dx2*sum_weight_b/delta)/(sum_weight_a+sum_weight_b*(delta+sum_weight_a)/delta)
    dbi = ((delta+sum_weight_a)*dai-dx2)/delta
    is_invalid = (x[i]-dxi <= b[j] && x[i] > b[j]) || (x[i] > a[j] && x[i]-dxi <= a[j]) || a[j]+dai > x_a_thd || x_a_thd == 1 || b[j]+dbi > x_b_thd || x_b_thd == 1
    return dai,dbi,is_invalid
end

function residual_update!(G::HyperGraphAndDegrees,
        x::Vector,a::Vector,b::Vector,i,dai,dbi,r,gamma,kappa,L::EdgeLoss,Q)
    Ht = G.Ht
    res_change = 0.0
    for k in Ht.colptr[i]:(Ht.colptr[i+1]-1)
        j = Ht.rowval[k]
        drij = loss_gradient(b[i]+dbi-x[j],L)-loss_gradient(x[j]-a[i]-dai,L)+loss_gradient(x[j]-a[i],L)-loss_gradient(b[i]-x[j],L)
        drij /= gamma
        rj = r[j]
        thd = kappa*G.deg[j]
        if rj <= thd && rj+drij > thd
            push!(Q,j)
        end
        r[j] += drij
        res_change += drij
    end
    return res_change
end

function _max_nz_degree(G::HyperGraphAndDegrees)
    H = G.H
    n = H.n
    maxd = zero(eltype(H.colptr))
    for i=1:n
        maxd = max(maxd, H.colptr[i+1]-H.colptr[i])
    end
    return maxd
end

function _max_edge_order(G::HyperGraphAndDegrees)
    Ht = G.Ht
    n = Ht.n
    maxo = zero(eltype(Ht.colptr))
    for i=1:n
        maxo = max(maxo, Ht.colptr[i+1]-Ht.colptr[i])
    end
    return maxo
end

"""
EdgeLoss{T} includes either TwoNormLoss or QHuberLoss, where we have
- `q` the value of q in the q-norm
- `delta` the value of delta in the q-Huber function
use loss_type(q,delta) for a type-unstable solution that will dispatch correctly

- `gamma` is for regularization, Infty returns seed set, 0 is hard/ill-posed.
- `kappa` is the sparsity regularilzation term.
- `rho` is the slack term in the KKT conditions to get faster convergence.
    (rho=1 is slow, rho=0)
- `eps` the value of epsilon in the local binary search
"""
function lh_diffusion(G::HyperGraphAndDegrees{T,Int64},S,gamma::T,kappa::T,rho::T,L::EdgeLoss{T};
        max_iters::Int=1000,x_eps::T=1.0e-8,aux_eps::T=1.0e-8,progress::Bool=true) where {T <: Real}
    init_x_eps,init_aux_eps = x_eps,aux_eps
    H,Ht = G.H,G.Ht
    n = size(Ht,1)
    m = size(Ht,2)
    x,a,b = zeros(n),zeros(m),zeros(m)
    """
    sum_weight is the sum of weight for xi
    a_thd is the smallest aj that is larger than xi
    b_thd is the smallest bj that is larger than xi
    sum_weight_a is the sum of weight for ai
    sum_weight_b is the sum of weight for bi
    x_a_thd is the smallest xj that is larger than ai
    x_b_thd is the smallest xj that is larger than bi
    """
    sums_weight,as_thd,bs_thd = copy(G.deg),ones(n),ones(n)
    sums_weight_a,sums_weight_b = Array{T}(copy(G.order)),Array{T}(copy(G.order))
    xs_a_thd,xs_b_thd = ones(m),ones(m)
    r = zeros(n)

    max_deg = _max_nz_degree(G)
    max_order = _max_edge_order(G)

    buf_x = zeros(max_order)
    buf_x_vals = zeros(max_order)
    buf_a = zeros(max_deg)
    buf_a_vals = zeros(max_deg)
    buf_b = zeros(max_deg)
    buf_b_vals = zeros(max_deg)
    Q = CircularDeque{Int}(n)
    for i in S
        r[i] = G.deg[i]
        push!(Q,i)
    end
    seedset = Set(S)

    iter = 0

    t0 = time()
    checkinterval = 10^4
    if progress == false
        checkinterval = max_iters
    end
    pushvol = 0
    nextcheck = checkinterval
    notify_time = 60*10.0
    last_time = t0
    last_iter = 0
    used_pm = false
    pm = Progress(max_iters, "LH: ")

    thd1 = minval(sum(G.deg[S])/sum(G.deg), L)
    thd2 = thd1
    # stop the loop early when the sum of residual continuously increase for this number of steps
    early_stop = 100
    stop = 0
    min_epsilon = min(x_eps,aux_eps)/10000 # the smallest episilon we will try

    while length(Q) > 0 && iter < max_iters && stop < early_stop
        min_dai_dbi = 1.0
        res_change = 0.0
        i = popfirst!(Q)
        if typeof(L) == TwoNormLoss{T}
            di = G.deg[i]
            dxi,ri_new,is_invalid = fast_dxi_solver!(x,r,i,di,rho,kappa,gamma,sums_weight,as_thd,bs_thd)
            # check if the assumption of fast_dxi_solver! holds
            if is_invalid
                dxi,ri_new,a_thd,b_thd,sum_weight = dxi_solver!(G,x,a,b,kappa,x_eps,gamma,r,seedset,rho,i,L,
                    buf_a,buf_b,buf_a_vals,buf_b_vals)
                as_thd[i],bs_thd[i] = a_thd,b_thd
                sums_weight[i] = sum_weight
            end
        else
            dxi,ri_new = dxi_solver!(G,x,a,b,kappa,x_eps,gamma,r,seedset,rho,i,L,
                buf_a,buf_b,buf_a_vals,buf_b_vals,thd1,thd2)
            thd2 = dxi/100.0
        end
        res_change += (ri_new - r[i])
        r[i] = ri_new
        x[i] += dxi
        for k in (H.colptr[i+1]-1):-1:H.colptr[i]
            j = H.rowval[k]
            if x[i]-dxi <= b[j] || x[i] > a[j]
                if typeof(L) != TwoNormLoss{T}
                    dai,dbi = dai_dbi_solver!(G,x,buf_x,buf_x_vals,a,b,aux_eps,j,L,dxi/(G.order[j]^(1/(L.q-1))))
                    # @show dai,dxi/(G.order[j]^(1/(L.q-1)))
                    # dai,dbi = dai_dbi_solver!(G,x,buf_x,buf_x_vals,a,b,aux_eps,j,L,dxi)
                else
                    dai,dbi,is_invalid = fast_dai_dbi_solver!(i,j,dxi,x,a,b,sums_weight_a,sums_weight_b,
                        xs_a_thd,xs_b_thd,G.delta)
                    # check if the assumption of fast_dai_dbi_solver! holds
                    if is_invalid
                        dai,dbi,x_a_thd,x_b_thd,sum_weight_a,sum_weight_b = dai_dbi_solver!(G,x,buf_x,buf_x_vals,a,b,aux_eps,j,L,dxi)
                        sums_weight_a[j],sums_weight_b[j] = sum_weight_a,sum_weight_b
                        xs_a_thd[j],xs_b_thd[j] = x_a_thd,x_b_thd
                    end
                end
                res_change += residual_update!(G,x,a,b,j,dai,dbi,r,gamma,kappa,L,Q)
                a[j] += dai
                b[j] += dbi
                min_dai_dbi = min(min_dai_dbi,min(dai,dbi))
            end
        end
        # adjust x_eps and aux_eps
        if dxi < x_eps*10 && x_eps > min_epsilon
            x_eps /= 10.0
        end
        if dxi > x_eps*100 && x_eps < init_x_eps
            x_eps *= 10.0
        end
        if min_dai_dbi < aux_eps*10 && aux_eps > min_epsilon
            aux_eps /= 10.0
        end
        if min_dai_dbi > aux_eps*100 && aux_eps < init_aux_eps && res_change < 0
            aux_eps *= 10.0
        end
        if min(x_eps,aux_eps) <= min_epsilon && res_change >= 0
            stop += 1
        else
            stop = 0
        end

        pushvol += H.colptr[i+1] - H.colptr[i]
        iter += 1

        if iter > nextcheck
            nextcheck = iter+checkinterval
            ct = time()

            if ct - t0 >= notify_time
                used_pm = true
                ProgressMeter.update!(pm, iter; showvalues =
                    [(:pushes_per_second,(iter-last_iter)/(ct-last_time)),
                     (:edges_per_second,pushvol/(ct-last_time))])
            end

            last_iter = iter
            last_time = ct
            pushvol = 0
        end
        # if iter % 1000 == 0
        #     # @show iter,i,sum(r),sum(x.>0),dxi,x_eps,aux_eps,stop
        #     @show iter,i,sum(r),dxi,thd1,thd2
        # end
        # @assert(res_change < 0)
        # @show iter
        # @show sparse(x).nzind,sparse(x).nzval
        # @show sparse(a).nzind,sparse(a).nzval
        # @show sparse(b).nzind,sparse(b).nzval
        # @show iter
        # @show r
        # @show Q,length(Q)
    end

    if used_pm == true
        ProgressMeter.finish!(pm)
    end

    if iter == max_iters && length(Q) > 0
        @warn "reached maximum iterations"
    end
    if stop >= early_stop && length(Q) > 0
        @warn "stopped early because sum of residual doesn't decrease for $early_stop steps"
    end
    return x,r,iter
end

end # end module

# include("common.jl")
# using Test
# @testset "LH" begin
#     name = "TinyZoo"
#     M = matread("small-hypergraphs/"*name*".mat")
#     H = M["H"]
#     order = vec(round.(Int64,sum(H,dims=2)))
#     good = findall(x->x > 0, order)
#     H = H[good,:]
#     G = LH.graph(H,1.0)
#     q = 1.8
#     delta = 0.0
#     L = LH.loss_type(q,delta)
#     kappa = 0.01
#     gamma = 0.1
#     rho = 0.999
#     S = [1]
#     x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
#     cond,cluster = hyper_sweepcut(G.H,x,G.deg,G.delta,0.0,G.order)
#     @test(abs(x[1]-0.1189) < 1.0e-4)
#     @test(abs(x[2]-0.01948) < 1.0e-4)
#     @test(abs(cond-0.12307692) < 1.0e-8)
#     true_cluster = [1, 4, 2, 18, 6, 5, 11, 10, 7]
#     @test(length(setdiff(Set(cluster),Set(true_cluster))) == 0)
#     q = 2.0
#     delta = 0.0
#     L = LH.loss_type(q,delta)
#     kappa = 0.01
#     gamma = 0.1
#     rho = 0.999
#     S = [1]
#     x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
#     @test(abs(x[1]-0.2009) < 1.0e-4)
#     @test(abs(x[2]-0.037891) < 1.0e-4)
#     G = LH.graph(H,100.0)
#     x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=10000)
#     @test(abs(x[1]-0.1338) < 1.0e-4)
#     @test(abs(x[2]-0.04130) < 1.0e-4)
#     q = 1.5
#     L = LH.loss_type(q,delta)
#     x,r,iter = LH.lh_diffusion(G,S,gamma,kappa,rho,L,max_iters=100000,x_eps=1.0e-16,aux_eps=1.0e-16)
#     @test(abs(x[1]-0.01184) < 1.0e-4)
#     @test(abs(x[2]-0.002087) < 1.0e-4)
# end
