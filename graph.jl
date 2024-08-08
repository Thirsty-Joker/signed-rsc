mutable struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int32, 1}
    v :: Array{Int32, 1} # uv is an edge
    w :: Array{Float64, 1} # weight of each edge
    #nbr :: Array{Array{Int32, 1}, 1}
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

# 有向图 度矩阵 
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

# 无向图 邻接矩阵
function adjacency_matrix2(G)
	u=zeros(2*G.m);
	v=zeros(2*G.m);
	w=zeros(2*G.m);
	for i=1:G.m
		u[i]=G.u[i];
		v[i]=G.v[i];
		w[i]=G.w[i];
		u[G.m+i]=G.v[i];
		v[G.m+i]=G.u[i];
		w[G.m+i]=G.w[i];
	end
	return sparse(u,v,w)
end

# 有向图 邻接矩阵
function adjacency_matrix_directed(G)
	return sparse(G.u,G.v,G.w,G.n,G.n)
end



function adjacency_matrix_pos(G :: Graph)
	F = spzeros(G.n, G.n)
	for i = 1 : G.m
		if G.w[i] > 0
			F[G.u[i], G.v[i]] += G.w[i]
			F[G.v[i], G.u[i]] += G.w[i]
		end
	end
	return F
end

function adjacency_matrix_pos2(G)
	u=zeros(2*G.m);
	v=zeros(2*G.m);
	w=zeros(2*G.m);
	for i=1:G.m
		if G.w[i]>0
			u[i]=G.u[i];
			v[i]=G.v[i];
			w[i]=G.w[i];
			u[G.m+i]=G.v[i];
			v[G.m+i]=G.u[i];
			w[G.m+i]=G.w[i];
		end
	end
	return sparse(u,v,w)
end

function adjacency_matrix_neg(G :: Graph)
	F = spzeros(G.n, G.n)
	for i = 1 : G.m
		if G.w[i] < 0
			F[G.u[i], G.v[i]] += G.w[i]
			F[G.v[i], G.u[i]] += G.w[i]
		end
	end
	return F
end

function adjacency_matrix_neg2(G)
	u=zeros(2*G.m);
	v=zeros(2*G.m);
	w=zeros(2*G.m);
	for i=1:G.m
		if G.w[i]<0
			u[i]=G.u[i];
			v[i]=G.v[i];
			w[i]=G.w[i];
			u[G.m+i]=G.v[i];
			v[G.m+i]=G.u[i];
			w[G.m+i]=G.w[i];
		end
	end
	return sparse(u,v,w)
end

function get_L(D,A)
	L=D-A
	return sparse(L)
end

#有向图 拉普拉斯矩阵
function laplacian_matrix(G :: Graph)
	d=zeros(G.n);
	for i=1:G.m
		# x=G.u[i];
		# y=G.v[i];
		# if G.w[i]>0
		# 	d[x]+=1;
		# 	d[y]+=1;
		# elseif G.w[i]<0
		# 	d[x]-=1;
		# 	d[y]-=1;
		# end
		d[G.u[i]]+=G.w[i];
	end
	# uu=zeros(2*G.m+G.n);
	# vv=zeros(2*G.m+G.n);
	# ww=zeros(2*G.m+G.n);
	uu=zeros(G.m+G.n);
	vv=zeros(G.m+G.n);
	ww=zeros(G.m+G.n);
	a=zeros(G.n);
	for i=1:G.n
		a[i]=i;
	end
	uu[1:G.m]=G.u;
	# uu[G.m+1:2*G.m]=G.v;
	uu[G.m+1:G.m+G.n]=a;
	vv[1:G.m]=G.v;
	# vv[G.m+1:2*G.m]=G.u;
	vv[G.m+1:G.m+G.n]=a;
	ww[1:G.m].=-G.w;
	# ww[G.m+1:2*G.m].=-G.w;
	ww[G.m+1:G.m+G.n]=d;
    return sparse(uu,vv,ww)
end

function get_2_H(D,A_pos,A_neg)
	H=[A_pos A_neg
	   A_neg A_pos]
	return sparse(H)
end

# 二倍图 邻接矩阵 无向
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

# 二倍图 邻接矩阵 有向
function twice_adjacency_directed(G)
	u=zeros(2*G.m)
	v=zeros(2*G.m)
	w=zeros(2*G.m)
	for i=1:G.m
		if G.w[i]>0
			u[i]=G.u[i]
			v[i]=G.v[i]
			w[i]=G.w[i]
			u[G.m+i]=G.n+G.u[i]
			v[G.m+i]=G.n+G.v[i]
			w[G.m+i]=G.w[i]
		elseif G.w[i]<0
			u[i]=G.u[i]
			v[i]=G.n+G.v[i]
			w[i]=-G.w[i]
			u[G.m+i]=G.n+G.u[i]
			v[G.m+i]=G.v[i]
			w[G.m+i]=-G.w[i]
		end
	end
	return sparse(u,v,w,2*G.n,2*G.n)
end

function get_2_D(D)
	n=size(D,1)
	D2 = zeros(2*n, 2*n)
    for i in 1:n
        for j in 1:n
            D2[i, j] = D[i, j]
            D2[i+n, j+n] = D[i, j]
        end
    end
    return D2
end

function get_2_L(D,A_pos,A_neg)
	L=[D-A_pos A_neg
	   A_neg D-A_pos]
	return sparse(L)
end


function get_B(G)
	F= spzeros(2*G.m,2*G.n)
	for i=1:G.m
		if G.w[i]>0
			F[i,G.u[i]]=1
			F[i,G.v[i]]=-1
			F[G.m+i,G.u[i]+G.n]=-1
			F[G.m+i,G.v[i]+G.n]=1
		elseif G.w[i]<0
			F[i,G.u[i]+G.n]=1
			F[i,G.v[i]]=-1
			F[G.m+i,G.u[i]]=1
			F[G.m+i,G.v[i]+G.n]=-1
		end
	end
	return F
end

# 关联矩阵
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

function getB_all(G)
	v1=zeros(2*G.m)
	v2=zeros(2*G.m)
	v1[1:G.m]=1:G.m;
	v1[G.m+1:2*G.m]=1:G.m;
	v2[1:G.m]=G.u;
	v2[G.m+1:2*G.m]=G.v;
	o=zeros(2*G.m)
	o[1:G.m].=1;
	o[G.m+1:2*G.m].=-G.w;
	return sparse(v1,v2,o)
end

function getB_pos(G)
	v1=zeros(2*G.m)
	v2=zeros(2*G.m)
	v1[1:G.m]=1:G.m;v1[G.m+1:2*G.m]=1:G.m;
	v2[1:G.m]=G.u;v2[G.m+1:2*G.m]=G.v;
	o=zeros(2*G.m)
	o[1:G.m].=1;o[G.m+1:2*G.m].=-max.(G.w,0);
	return sparse(v1,v2,o)
end

function getB_neg(G)
	v1=zeros(2*G.m)
	v2=zeros(2*G.m)
	v1[1:G.m]=1:G.m;v1[G.m+1:2*G.m]=1:G.m;
	v2[1:G.m]=G.u;v2[G.m+1:2*G.m]=G.v;
	o=zeros(2*G.m)
	o[1:G.m].=1;o[G.m+1:2*G.m].=-min.(G.w,0);
	return sparse(v1,v2,o)
end


function Uniform(n)
    x = rand(n)
    return x
end