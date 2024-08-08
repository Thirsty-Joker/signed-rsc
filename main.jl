include("new_graph.jl")
include("algorithm.jl")

using GraphDatasets
using GeneralGraphs, LinearAlgebraUtils
using LinearAlgebra, SparseArrays

function initialize(file::AbstractString)
    file=split(file)[1]
    g_signed=NormalWeightedGraph{Float64}(file,"data/"*file)
    _,A_signed=diagadj(g_signed)
    degree,alpha=create_new_graph(file)
    println("create new graph done")
    newfile="new-data/"*file[1:end-4]*"-new.txt"
    g = NormalWeightedDiGraph{Float64}(file, newfile)
    n= round(Int,num_nodes(g)/2)
    println("n=",n)
    println("alpha=",alpha)
    return A_signed,g,degree,alpha,n
end

fname=open("input/data.txt","r")
str=readline(fname)
num=parse(Int,str)
for nnnn in 1:num
    println("Number=$nnnn")
    str=readline(fname)
    A_signed,g,degree,alpha,n=initialize(str)
    println("start computing")

    #exact
    q1=zeros(n)
    std=inv(I - alpha* Matrix(A_signed))
    for o in 1:n
        q1[o]=std[o,o]
    end

    G=read_data(file)
    q2=zeros(n)
    q3=zeros(n)

    eps=[0.3,0.2,0.1]
    lll=[2000,6000,10000]

    for m in 1:1
        ll=lll[m]
        epsilon=eps[m]

        #approx_forest
        c2 = simple_wilson(g, ll) 
        q2=subgraph_centrality(c2,g,alpha,degree)

        #enhanced_forest
        c3=simple_enhanced_wilson(g,ll)
        q3=subgraph_centrality(c3,g,alpha,degree)

        #solver
        q4=jl_solver(alpha,G,degree,epsilon)

        for i in 1:n
            println(q1[i]," ",q2[i]," ",q3[i]," ",q4[i])
        end
    end    
    println("---------------------------------------------")
end

