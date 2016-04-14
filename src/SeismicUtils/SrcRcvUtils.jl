
export readSrcRcvLocationFile,generateSrcRcvProjOperators,readDataFileToDataMat,writeSrcRcvLocFile,writeDataFile,limitDataToOffset

function readSrcRcvLocationFile(filename,M::RegularMesh)
# read files for sources and receivers location in UTM and translate them into indices locations.
# file contains a table: idx xUTM yUTM zUTM
A = readdlm(filename);
domainBoundaryUTM = M.domain;
n = M.n;
# domainBoundaryUTM: [x_start, x_end, y_start, y_end, z_start, z_end];
# n: number of cells. (n-1: number of cells)

for k=1:M.dim
	A[:,k+1] = A[:,k+1] - domainBoundaryUTM[2*k-1];
	A[:,k+1] = A[:,k+1] ./ (domainBoundaryUTM[2*k] - domainBoundaryUTM[2*k-1]);
	A[:,k+1] = round(A[:,k+1] .* (n[k]))+1.0;
end
A = round(Int,A);
return A;
end

function generateSrcRcvProjOperators(n_nodes::Vector,SrcRcvNodeLoc::Array{Int64})
# gets a sources/receivers location table (in indices) and generates operators.
nSrcRcv = size(SrcRcvNodeLoc,1);
I = zeros(Int64,nSrcRcv);
J = zeros(Int64,nSrcRcv);
V = zeros(Float64,nSrcRcv);
for k=1:nSrcRcv
	J[k] = k;
	if length(n_nodes)==2
		I[k] = sub2ind(tuple(n_nodes...),SrcRcvNodeLoc[k,2],SrcRcvNodeLoc[k,3]);
	else
		I[k] = sub2ind(tuple(n_nodes...),SrcRcvNodeLoc[k,2],SrcRcvNodeLoc[k,3],SrcRcvNodeLoc[k,4]);
	end
	V[k] = 1.0;
end
P = sparse(I,J,V,prod(n_nodes),nSrcRcv);
end

function readDataFileToDataMat(filename,SrcNodeLoc::Array{Int64},RcvNodeLoc::Array{Int64})
# reads a data file into a data matrix D for the code to work with 
Data = readdlm(filename);
# Data is a table of #src #rcv tau(src,rcv) Wd(src,rcv) and the rest is zero.
# CellLocMap is a table of #src i j k 
numSrc = size(SrcNodeLoc,1);
numRcv = size(RcvNodeLoc,1);

if size(Data,2) == 4
	D = zeros(Float64,numRcv,numSrc);
	Wd = zeros(Float64,numRcv,numSrc);
	SrcIdxs = SrcNodeLoc[:,1];
	RcvIdxs = RcvNodeLoc[:,1];
	for k=1:size(Data,1)
		srcIdx = convert(Int64,round(Data[k,1]));
		rcvIdx = convert(Int64,round(Data[k,2]));
		isrc = binarySearch(SrcIdxs,srcIdx);
		ircv = binarySearch(RcvIdxs,rcvIdx);
		D[ircv,isrc] = Data[k,3] ;#+ 1im*Data[k,4]; 
		Wd[ircv,isrc] = Data[k,4] ;#+ 1im*Data[k,5];
	end
elseif size(Data,2) == 6
	D = zeros(Complex128,numRcv,numSrc);
	Wd = zeros(Complex128,numRcv,numSrc);
	SrcIdxs = SrcNodeLoc[:,1];
	RcvIdxs = RcvNodeLoc[:,1];
	for k=1:size(Data,1)
		srcIdx = convert(Int64,round(Data[k,1]));
		rcvIdx = convert(Int64,round(Data[k,2]));
		isrc = binarySearch(SrcIdxs,srcIdx);
		ircv = binarySearch(RcvIdxs,rcvIdx);
		D[ircv,isrc] = Data[k,3] + 1im*Data[k,5]; 
		Wd[ircv,isrc] = Data[k,4] + 1im*Data[k,6];
	end

end
return D,Wd;
end

function binarySearch(arr, value)
    low = 1
    high = length(arr)
    while low <= high
        mid = round(Int,(low+high)/2)
        if arr[mid] > value 
          high = mid-1
        elseif arr[mid] < value 
          low = mid+1
        else 
          return mid
        end
    end     
    return -1
end

function writeSrcRcvLocFile(filename,Msh::RegularMesh,pad::Int64,jump::Int64)
# writes an equally distanced sources/receivers file.
n_tup = tuple((Msh.n+1)...); # n_tup is nodes...
Q = zeros(n_tup);

if Msh.dim==2
	Q[pad+1:jump:end-pad-1,1] = 1.0
else	
	Q[pad+1:jump:end-pad-1,pad+1:jump:end-pad-1,1] = 1.0;
end
numSrcRcv = round(Int,sum(Q));
Q = find(Q);
SrcRcvTable = zeros(numSrcRcv,1+Msh.dim);
SrcRcvIDs = randperm(numSrcRcv);
SrcRcvIDs = sort(SrcRcvIDs[1:numSrcRcv]);
SrcRcvTable[:,1] = SrcRcvIDs;
startLocations = Msh.domain[1:2:end];
for k=1:numSrcRcv
	SrcRcvTable[k,2:end] = (collect(ind2sub(n_tup,Q[k]))-1).*Msh.h + startLocations;
end
writedlm(filename,SrcRcvTable);
end


function writeDataFile(filename,D::Union{Array{Float64,2},Array{Complex128,2}},Wd::Union{Array{Float64,2},Array{Complex128,2}},SrcNodeLoc::Array{Int64},RcvNodeLoc::Array{Int64})
srcIDs = SrcNodeLoc[:,1];
rcvIDs = RcvNodeLoc[:,1];

if eltype(D)==Float64
	println("Writing real data file");  
	Data = zeros(Float32,countnz(Wd),4);
	k = 1;
	for j = 1:size(D,2) 
		for i = 1:size(D,1)
			if Wd[i,j] != 0.0
				Data[k,1] = srcIDs[j];
				Data[k,2] = rcvIDs[i];
				Data[k,3] = D[i,j];
				Data[k,4] = Wd[i,j];
				k+=1;
			end
		end
	end
else
	println("Writing ``complex'' data file");  
	Data = zeros(Float32,countnz(Wd),6);
	k = 1;
	for j = 1:size(D,2) 
		for i = 1:size(D,1)
			if Wd[i,j] != 0.0
				Data[k,1] = srcIDs[j];
				Data[k,2] = rcvIDs[i];
				Data[k,3] = real(D[i,j]);
				Data[k,4] = real(Wd[i,j]);
				Data[k,5] = imag(D[i,j]);
				Data[k,6] = imag(Wd[i,j]);
				k+=1;
			end
		end
	end
end
writedlm(filename,Data);
return;
end

function limitDataToOffset(Wd::Union{Array{Float64,2},Array{Complex128,2}},SrcNodeLoc::Array{Int64},RcvNodeLoc::Array{Int64},offsetInNodes::Int64)
for j = 1:size(Wd,2)
	srcLoc = SrcNodeLoc[j,2:end];
	for i = 1:size(Wd,1)
		rcvLoc = RcvNodeLoc[i,2:end];
		dist = norm(srcLoc - rcvLoc);
		# dist = rcvLoc[1] - srcLoc[1];
		if dist > offsetInNodes || dist < 2
			Wd[i,j] = 0.0;
		end
	end
end
return Wd;
end

