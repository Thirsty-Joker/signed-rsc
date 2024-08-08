include("graph.jl")

function create_new_graph(filename::AbstractString)
    G=read_data(filename)
    println("read data done")

    println("------",G.m,"------")

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

    u=zeros(Int,2*G.m)
    v=zeros(Int,2*G.m)
    w=zeros(Float32,2*G.m)

    for i=1:G.m
        degree[G.u[i]]+=G.w[i]
        degree[G.v[i]]+=G.w[i]
    end

    maxdegree=0
    for i in 1:G.n
        if degree[i]>maxdegree
            maxdegree=degree[i]
        end
    end
    alpha=Float32(1.0/(maxdegree+1))

    for i=1:G.m
        u[i]=G.u[i]
        v[i]=G.v[i]
        w[i]=alpha/(1-alpha*degree[G.u[i]])
        u[i+G.m]=G.v[i]
        v[i+G.m]=G.u[i]
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