if nworkers()<=1
	addprocs(2);
end

using  jInv.Mesh
using  jInv.Utils
using  jInv.InverseSolve
using  EikonalInv
using  MAT

#############################################################################################################
plotting = false;
if plotting
	using  PyPlot
	close("all");
end
#############################################################################################################
dataDir = pwd();
resultsDir = pwd();

###############################################################################################################
# this filename can be used to save the model during the iterations. If set to "", no images or files will be saved.
resultsFilenamePrefix = "";
# resultsFilenamePrefix = string(resultsDir,"/travelTimeInvSEG"); 
################################################################################################################


modelDir = pwd();
include("../drivers/prepareTravelTimeDataFiles.jl");
include("../drivers/plotModel.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");
include("../drivers/setupTravelTimeTomography.jl");

########################################################################################################
######################################## for SEG 256X512 ###############################################
#######################################################################################################
dim     = 2;
pad     = 0;
jump    = 10;
offset  = 256;
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[256,128],[0.0,13.5,0.0,4.2]);
dataFilenamePrefix = string(dataDir,"/DATA_SEG",tuple((Minv.n+1)...));


if plotting
	limits = (minimum(m),maximum(m));
	figure()
	plotModel(m,true,true,Minv,pad,limits);
	figure()
	plotModel(mref,true,true,Minv,pad,limits);
end

prepareTravelTimeDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,pad,jump,offset);


(Q,P,pMisRFs,SourcesSubInd,contDiv,Iact,mback,mref,boundsHigh,boundsLow,resultsFilename) = setupTravelTimeTomography(m,dataFilenamePrefix, resultsFilenamePrefix);


########################################################################################################
##### Set up Inversion #################################################################################
########################################################################################################

maxStep=0.2*maximum(boundsHigh);

a = minimum(boundsLow)*0.8;
b = maximum(boundsHigh)*1.2;
modfun(x) = getBoundModel(x,a,b);
mref = getBoundModelInv(mref,a,b);
boundsHigh = boundsHigh*10000000.0;
boundsLow = -boundsLow*100000000.0

cgit  = 10;
maxit = 7;
alpha = 1e+0;
pcgTol = 1e-1;

HesPrec=getSSORCGRegularizationPreconditioner(1.0,1e-5,1000)

################################################# GIT VERSION OF JINV #################################################

regparams = [1.0,1.0,1.0,1e-2];
regfun(m, mref, M) = wdiffusionRegNodal(m, mref, M, Iact=Iact, C = regparams);
	
	
function dump(mc,Dc,iter,pInv,pMis)
	if resultsFilename!=""
		fullMc = reshape(Iact*modfun(mc)[1] + mback,tuple((pInv.MInv.n+1)...));
		Temp = splitext(resultsFilename);
		Temp = string(Temp[1],"_GN",iter,Temp[2]);
		writedlm(Temp,convert(Array{Float16},fullMc));
		if plotting
			close(888);
			figure(888);
			println(splitdir(Temp)[2])
			plotModel(fullMc,true,false,[],0,[a,b],splitdir(Temp)[2]);
		end
	end
end	

pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
					 
mc,Dc,flag = projGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);

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

rm("DATA_SEG(256,128)_travelTime.dat");
rm("DATA_SEG(256,128)_rcvMap.dat");
rm("DATA_SEG(256,128)_srcMap.dat");
rm("DATA_SEG(256,128)_PARAM.mat");
rm("jInv.out");
#############################################################################################

