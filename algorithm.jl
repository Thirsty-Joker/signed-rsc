using GraphDatasets
using GeneralGraphs, LinearAlgebraUtils, ProgressMeter
using LinearAlgebra, SparseArrays
using StatsBase, Random

function simple_wilson(g::NormalWeightedDiGraph{Float64}, sample_num::Int)
    n = num_nodes(g)
    node=round(Int,n/2)
    in_forests = Vector{Bool}(undef, n)
    next = Vector{Int}(undef, n)
    root = Vector{Int}(undef, n)
    wei=zeros(Float64,n)
    for u in 1:n
        wei[u]=sum(g.weights[u])
    end
    ans = zeros(Float64, n)
    for _ in 1:sample_num
        fill!(in_forests, false)
        for src in 1:node
            u = src
            while in_forests[u] == false
                if rand(Float64) * (wei[u] + 1) < 1
                    in_forests[u] = true
                    root[u] = u
                    ans[u] += 1/sample_num
                    break
                end
                next[u] = StatsBase.sample(g.adjs[u], Weights(g.weights[u]))
                u = next[u]
            end
            r = root[u]
            u = src
            while in_forests[u] == false
                in_forests[u] = true
                root[u] = r
                if r==u+node
                    ans[u]-=1/sample_num
                end
                u = next[u]
            end
        end
    end
    return ans
end

function simple_enhanced_wilson(g::NormalWeightedDiGraph{Float64}, sample_num::Int)
    n = num_nodes(g)
    nodes= round(Int,n/2)
    _,sp_A=diagadj(g)
    in_forests = Vector{Bool}(undef,n)
    next = Vector{Int}(undef,n)
    root = Vector{Int}(undef,n)
    wei=zeros(Float64,n)
    ans = zeros(Float64, n)
    for u in 1:n
        wei[u]=sum(g.weights[u])
        ans[u]=1/(1+wei[u])
    end
    for _ in 1:sample_num
        fill!(in_forests, false)
        for src in 1:nodes
            u = src
            while in_forests[u] == false
                if rand(Float64) * (wei[u] + 1) < 1
                    in_forests[u] = true
                    root[u] = u
                    break
                end
                next[u] = StatsBase.sample(g.adjs[u], Weights(g.weights[u]))
                u = next[u]
            end
            r = root[u]
            u = src
            while in_forests[u] == false
                in_forests[u] = true
                root[u] = r
                if u<=nodes
                    ru=sp_A[r,u]
                    runode=sp_A[r,u+nodes]
                    if ru!=0
                        ans[u]+=ru/((1+wei[u])*sample_num)
                    end
                    if runode!=0
                        ans[u]-=runode/((1+wei[u])*sample_num)
                    end
                end
                u = next[u]
            end
        end
    end
    return ans
end

function jl_solver(alpha,G,degree,epsilon)
    n=G.n
    m=G.m
    A2=twice_adjacency(G)
    B=incidence_matrix(G)
    u=zeros(2*n)
    for i in 1:2*n
        u[i]=i
    end
    q=zeros(2*n)
    for i in 1:2*n
        q[i]=1-alpha*degree[i]
    end
    Q=sparse(u,u,q)
    k = round(Int64, 24*log(G.n) / (epsilon^2)) + 1
    println("k=$k")
    c=zeros(n)
    for _ in 1:k
        m1=rand([1,-1],2*n)
        p1=rand([1,-1],2*m)
        f=approxchol_sddm(I-alpha*A2)
        x=sqrt.(Q)*m1
        y=B'*p1
        z1=f(x)
        z2=f(y)
        for i in 1:n
            c[i]+=((z1[i]-z1[i+n])^2/k+alpha*(z2[i]-z2[i+n])^2/k)/2
        end
    end
    return c
end


function subgraph_centrality(c::Vector{Float64}, g::NormalWeightedDiGraph{Float64}, a::Float32, degree::Vector{Float32})
    n = num_nodes(g)
    q = zeros(Float64,n)
    for i=1:round(Int,n/2)
        q[i] = c[i]/ (1 - a*degree[i])
    end
    return q[1:round(Int,n/2)]
end
