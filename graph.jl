mutable struct Graph
    n :: Int 
    m :: Int 
    u :: Array{Int32, 1}
    v :: Array{Int32, 1} 
    w :: Array{Float64, 1} 
end
using LinearAlgebra
using SparseArrays
using Laplacians


function read_data(filename::AbstractString)
    ffname=string("data/",filename)
    open(ffname) do file
        G=Graph(0,0,Int[],Int[],Float32[])
        for line in eachline(file)
            fields = split(line)
            if parse(Float32, fields[3])>0
                push!(G.u, parse(Int, fields[1]))
                push!(G.v, parse(Int, fields[2]))
                push!(G.w, parse(Float32, fields[3]))
                G.m+=1
            elseif parse(Float32, fields[3])<=0
                push!(G.u, parse(Int, fields[1]))
                push!(G.v, parse(Int, fields[2]))
                push!(G.w, parse(Float32, fields[3]))
                G.m+=1
            end
            if eof(file)
                break
            end
        end
        G.n=max(maximum(G.u),maximum(G.v))
        return G
    end
end

function degree_matrix2(G)
	u=zeros(G.n);
	d=zeros(G.n);
	for i=1:G.n
		u[i]=i
	end
	for i=1:G.m
		d[G.u[i]]+=G.w[i];
	end
	return sparse(u,u,d)
end

function twice_adjacency(G)
	u=zeros(4*G.m)
	v=zeros(4*G.m)
	w=zeros(4*G.m)
	for i=1:G.m
		if G.w[i]>0
			u[i]=G.u[i]
			v[i]=G.v[i]
			w[i]=G.w[i]
			u[G.m+i]=G.n+G.u[i]
			v[G.m+i]=G.n+G.v[i]
			w[G.m+i]=G.w[i]
			u[2*G.m+i]=G.v[i]
			v[2*G.m+i]=G.u[i]
			w[2*G.m+i]=G.w[i]
			u[3*G.m+i]=G.n+G.v[i]
			v[3*G.m+i]=G.n+G.u[i]
			w[3*G.m+i]=G.w[i]
		elseif G.w[i]<0
			u[i]=G.u[i]
			v[i]=G.n+G.v[i]
			w[i]=-G.w[i]
			u[G.m+i]=G.n+G.u[i]
			v[G.m+i]=G.v[i]
			w[G.m+i]=-G.w[i]
			u[2*G.m+i]=G.v[i]
			v[2*G.m+i]=G.n+G.u[i]
			w[2*G.m+i]=-G.w[i]
			u[3*G.m+i]=G.n+G.v[i]
			v[3*G.m+i]=G.u[i]
			w[3*G.m+i]=-G.w[i]
		end
	end
	return sparse(u,v,w)
end

function incidence_matrix(G)
	u=zeros(4*G.m)
	v=zeros(4*G.m)
	w=zeros(4*G.m)
	for i=1:G.m
		if G.w[i]>0
			u[i]=i
			v[i]=G.u[i]
			w[i]=1
			u[G.m+i]=i
			v[G.m+i]=G.v[i]
			w[G.m+i]=-1
			u[2*G.m+i]=i+G.m
			v[2*G.m+i]=G.n+G.u[i]
			w[2*G.m+i]=-1
			u[3*G.m+i]=i+G.m
			v[3*G.m+i]=G.n+G.v[i]
			w[3*G.m+i]=1
		elseif G.w[i]<0
			u[i]=i
			v[i]=G.u[i]+G.n
			w[i]=1
			u[G.m+i]=i
			v[G.m+i]=G.v[i]
			w[G.m+i]=-1
			u[2*G.m+i]=i+G.m
			v[2*G.m+i]=G.u[i]
			w[2*G.m+i]=1
			u[3*G.m+i]=i+G.m
			v[3*G.m+i]=G.v[i]+G.n
			w[3*G.m+i]=-1
		end
	end
	return sparse(u,v,w)
end