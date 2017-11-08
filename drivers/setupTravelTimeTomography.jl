function setupTravelTimeTomography(m,filenamePrefix::String, resultsOutputFolderAndPrefix::String,useFilesForFields::Bool = false)


file = matopen(string(filenamePrefix,"_PARAM.mat"));
n_cells = read(file,"MinvN");
OmegaDomain = read(file,"MinvOmega");
Minv = getRegularMesh(OmegaDomain,n_cells);
HO = read(file,"HO");

boundsLow = read(file,"boundsLow");
boundsHigh = read(file,"boundsHigh");
mref =  read(file,"mref");
close(file);

resultsFilename = "";
if resultsOutputFolderAndPrefix!=""
	resultsFilename = string(resultsOutputFolderAndPrefix,tuple((Minv.n+1)...),".dat");
end
###########################################################

### Read receivers and sources files
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");

srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);

println("Travel time tomography: ",size(Q,2)," sources.");
#############################################################################################


println("Reading data:");

(DobsEik,WdEik) =  readDataFileToDataMat(string(filenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap);

N = prod(Minv.n+1);

Iact = speye(Float16,N);
Iact = convert(SparseMatrixCSC{Float16,Int32},Iact);
mback   = zeros(Float64,N);


## Setting the sea constant:
mask = zeros(N);
sea = abs.(m[:] .- minimum(m)) .< 5e-2;
mask[sea] = 1.0;
# # setup active cells
mback = vec(m[:].*mask);
sback = velocityToSlowSquared(mback)[1];
sback[mask .== 0.0] = 0.0;

Iact = Iact[:,mask .== 0.0];
boundsLow = Iact'*boundsLow;
boundsHigh = Iact'*boundsHigh;
mref = Iact'*mref[:];


########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

EikMPIWorkers = nworkers(); # this just set the maximal MPI workers. To activate parallelism, run addprocs()

(pFor,contDiv,SourcesSubInd) = getEikonalInvParam(Minv,Q,P,HO,EikMPIWorkers,useFilesForFields);


misfun = SSDFun

Wd 		= Array{Array{Float64}}(length(pFor));
dobs 	= Array{Array{Float64}}(length(pFor));
for i=1:length(pFor)
	I_i = SourcesSubInd[i];
	Wd[i]   = WdEik[:,I_i];
	dobs[i] = DobsEik[:,I_i];
end
WdEik = 0;
# DobsEik = 0;

pMisRFs = getMisfitParam(pFor, Wd, dobs, misfun, Iact,sback);

return (Q,P,pMisRFs,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename,DobsEik);
end


