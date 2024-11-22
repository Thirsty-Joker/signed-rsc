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
for file_num in 1:num
    @show file_num
    str=readline(fname)
    A_signed,g,degree,alpha,n=initialize(str)
    println("start computing")

    #exact
    exactRSC=zeros(n)
    std=inv(I - alpha* Matrix(A_signed))
    for o in 1:n
        exactRSC[o]=std[o,o]
    end

    G=read_data(file)

    eps=[0.2,0.1]
    lll=[1000,2000]

    for m in 1:2
        ll=lll[m]
        epsilon=eps[m]

        #simpleRSC
        simpleRSC = simple_RSC(g, ll) 

        #groupRSC
        groupRSC=group_RSC(g,ll)

        #neighborRSC
        neighborRSC=neighbor_RSC(g,ll)

        #solverRSC
        solverRSC=jl_solver(alpha,G,degree,epsilon)

        @show simpleRSC,groupRSC,neighborRSC,solverRSC
    end    
    println("---------------------------------------------")
end

