
export arrangeRemoteCallDataIntoLocalData,getSourcesIndicesOfKthWorker

function arrangeRemoteCallDataIntoLocalData(DremoteCall::Union{Array{Array{Float64,2},1},Array{Array{Complex128,2},1}})
N = length(DremoteCall);
l = 0;
for k=1:N
	l+=size(DremoteCall[k],2);
end
Dlocal = zeros(eltype(DremoteCall[1]),size(DremoteCall[1],1),l)
for k=1:N
	Dlocal[:,k:N:end] = DremoteCall[k];
end	
return Dlocal;
end


function getSourcesIndicesOfKthWorker(nworkers::Int64,k::Int64,nsrc::Int64)
return k:nworkers:nsrc;
end
