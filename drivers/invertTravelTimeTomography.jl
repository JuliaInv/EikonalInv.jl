function invertTravelTimeTomography(m,filenamePrefix::ASCIIString, resultsOutputFolderAndPrefix::ASCIIString,maxit = 10)

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
sea = abs(m[:] .- maximum(m)) .< 1e-2;
mask[sea] = 1;
# # setup active cells
mback = vec(m[:].*mask);
Iact = Iact[:,mask .== 0.0];
boundsLow = Iact'*boundsLow;
boundsHigh = Iact'*boundsHigh;
mref = Iact'*mref[:];


########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

EikMPIWorkers = nworkers(); # this just set the maximal MPI workers. To activate parallelism, run addprocs()

(pFor,contDiv,SourcesSubInd) = getEikonalInvParam(Minv,Q,P,HO,EikMPIWorkers);


misfun = SSDFun

Wd 		= Array(Array{Float64},length(pFor));
dobs 	= Array(Array{Float64},length(pFor));
for i=1:length(pFor)
	I_i = SourcesSubInd[i];
	Wd[i]   = WdEik[:,I_i];
	dobs[i] = DobsEik[:,I_i];
end
WdEik = 0;
DobsEik = 0;

pMisRFs = getMisfitParam(pFor, Wd, dobs, misfun, Iact,mback);


########################################################################################################
##### Set up Inversion #################################################################################
########################################################################################################


maxStep=0.2*maximum(boundsHigh);

a = minimum(boundsLow)*0.8;
b = maximum(boundsHigh)*1.2;

modfun = identityMod;


cgit = 10; 
alpha = 1e+0;
pcgTol = 1e-1;

HesPrec=getSSORCGRegularizationPreconditioner(1.0,1e-5,1000)

################################################# GIT VERSION OF JINV #################################################

regparams = [1.0,1.0,1.0,1e-2];
regfun(m, mref, M) = wdiffusionRegNodal(m, mref, M, Iact=Iact, C = regparams);
	
pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
function dump(mc,Dc,iter,pInv,pMis)
	if resultsFilename!=""
		fullMc = reshape(Iact*modfun(mc)[1] + mback,tuple((pInv.MInv.n+1)...));
		Temp = splitext(resultsFilename);
		Temp = string(Temp[1],"_GN",iter,Temp[2]);
		writedlm(Temp,convert(Array{Float16},fullMc));
		if plotting
			close(888);
			figure(888);
			plotModel(fullMc,true,false,[],0,[a,b],splitdir(Temp)[2]);
		end
	end
end						 
						 
tic()
mc,Dc,flag = projGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);
# mc,Dc,flag = ADMMconst(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);
toc()

Dpred = Array(Array{Float64,2},length(pMisRFs))
for k = 1:length(pMisRFs)
	Dpred[k] = fetch(Dc[k]);
end
Dpred = arrangeRemoteCallDataIntoLocalData(Dpred);

if resultsFilename!=""
	Temp = splitext(resultsFilename);
	writedlm(string(Temp[1],"_predictedData",Temp[2]),convert(Array{Float16},Dpred));
	writedlm(string(Temp[1],"_recoveredModel",Temp[2]),convert(Array{Float16},reshape(Iact*modfun(mc)[1] + mback,tuple((pInv.MInv.n+1)...))));
end
########################################################################################################################
return mc,Dpred;
end