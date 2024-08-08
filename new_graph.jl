include("graph.jl")
# using GraphDatasets
# using GeneralGraphs, LinearAlgebraUtils, ProgressMeter
# using LinearAlgebra, SparseArrays

#无向图 ---> 有向有权图
function create_new_graph(filename::AbstractString)
    G=read_data(filename)
    println("read data done")

    println("------",G.m,"------")

    #符号图->二倍无符号图
    u=zeros(Int,2*G.m)
    v=zeros(Int,2*G.m)
    w=zeros(Float32,2*G.m)
    degree=zeros(Float32,2*G.n)
    for i=1:G.m   
        if G.w[i]>=0
            u[i]=G.u[i]
            v[i]=G.v[i]
            w[i]=G.w[i]
            u[i+G.m]=G.u[i]+G.n
            v[i+G.m]=G.v[i]+G.n
            w[i+G.m]=G.w[i]
        elseif G.w[i]<0
            u[i]=G.u[i]
            v[i]=G.v[i]+G.n
            w[i]=-G.w[i]
            u[i+G.m]=G.u[i]+G.n
            v[i+G.m]=G.v[i]
            w[i+G.m]=-G.w[i]
        end
    end
    G.n=2*G.n
    G.m=2*G.m
    G.u=u
    G.v=v
    G.w=w
    println("unsigned done")

    #无向图->有向图
    u=zeros(Int,2*G.m)
    v=zeros(Int,2*G.m)
    w=zeros(Float32,2*G.m)


    for i=1:G.m
        degree[G.u[i]]+=G.w[i]
        degree[G.v[i]]+=G.w[i]
    end

    # D=zeros(Float32,G.n,G.n)
    # for i=1:G.n
    #     D[i,i]=degree[i]
    # end

    maxdegree=0
    for i in 1:G.n
        # degree[i]=D[i,i]
        if degree[i]>maxdegree
            maxdegree=degree[i]
        end
    end
    # alpha=Float32(1/maxdegree-0.00001)
    alpha=Float32(1.0/(maxdegree+1))
    println("最大度max_degree=",maxdegree)
    # alpha=Float32(0.25)

    for i=1:G.m
        u[i]=G.u[i]
        v[i]=G.v[i]
        # w[i]=alpha/(1-alpha*D[G.u[i],G.u[i]])
        w[i]=alpha/(1-alpha*degree[G.u[i]])
        u[i+G.m]=G.v[i]
        v[i+G.m]=G.u[i]
        # w[i+G.m]=alpha/(1-alpha*D[G.v[i],G.v[i]])
        w[i+G.m]=alpha/(1-alpha*degree[G.v[i]])
    end
    println("new graph done")

    file = open("new-data/"*filename[1:end-4]*"-new.txt", "w")

    println("-----",length(u),"------")

    for i in 1:length(u)
        write(file, "$(u[i]) $(v[i]) $(w[i])\n")
    end
    close(file)
    println("write file done")


    return degree,alpha
end