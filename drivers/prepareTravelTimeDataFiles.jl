function prepareTravelTimeDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,filenamePrefix::ASCIIString,pad::Int64,jump::Int64,offset::Int64)
### Here we generate data files for sources and receivers. This can be replaced with real sources/receivers files.
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,pad,jump);
writeSrcRcvLocFile(RCVfile,Minv,pad,1);

dataFullFilename = string(filenamePrefix,"_travelTime");
HO = false;
prepareTravelTimeDataFiles(m, Minv, filenamePrefix,dataFullFilename,offset,HO);


file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"MinvOmega",Minv.domain);
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"MinvN",Minv.n);
write(file,"HO",HO);
close(file);
return;
end

function prepareTravelTimeDataFiles(m, Minv::RegularMesh, filenamePrefix::ASCIIString,dataFullFilename::ASCIIString, offset::Int64,HO::Bool)

RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*(1/(norm(Minv.h)^2));
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);

# compute observed data

println("~~~~~~~ Getting data Eikonal: ~~~~~~~");
(pForEIK,contDivEIK,SourcesSubIndEIK) = getEikonalInvParam(Minv,Q,P,HO,nworkers());

Mesh2Mesh = Array(RemoteRef{Channel{Any}},length(pForEIK))	
for i=1:length(Mesh2Mesh)
	k = pForEIK[i].where
	Mesh2Mesh[i] = remotecall(k,speye,prod(Minv.n+1))
end

(D,pForEIK) = getData(m[:],pForEIK,Mesh2Mesh,true);


Dobs = Array(Array{Float64,2},length(pForEIK))
for k = 1:length(pForEIK)
	Dobs[k] = fetch(D[k]);
end
Dobs = arrangeRemoteCallDataIntoLocalData(Dobs);

# D should be of length 1 becasue constMUSTBeOne = 1;
# Dobs = fetch(D[1]);
Dobs += 0.01*mean(abs(Dobs))*randn(size(Dobs,1),size(Dobs,2));
Wd = (50.0./(abs(Dobs)+ 0.1*mean(abs(Dobs))));
Wd = limitDataToOffset(Wd,srcNodeMap,rcvNodeMap,offset);
writeDataFile(string(dataFullFilename,".dat"),Dobs,Wd,srcNodeMap,rcvNodeMap);

return (Dobs,Wd)
end