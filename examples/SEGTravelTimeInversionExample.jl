if nworkers()<=1
	println("Make code faster by adding processes: addprocs(X)")
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
jump    = 5;
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

(Q,P,pMisRFs,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename) = setupTravelTimeTomography(m,dataFilenamePrefix, resultsFilenamePrefix);

########################################################################################################
##### Set up Inversion #################################################################################
########################################################################################################

maxStep=0.1*maximum(boundsHigh);
modfun = velocityToSlowSquared;
a = minimum(boundsLow);
b = maximum(boundsHigh);

################################################################################################################
############### USE A BOUND MODEL ##############################################################################
################################################################################################################
# maxStep=0.1*maximum(boundsHigh);
# a = minimum(boundsLow);
# b = maximum(boundsHigh);
# function modfun(m)
	# bm,dbm = getBoundModel(m,a,b);
	# bs,dbs = velocityToSlowSquared(bm);
	# dsdm = dbm*dbs;
	# return bs,dsdm;
# end
# mref = getBoundModelInv(mref,a,b);
# boundsHigh = boundsHigh*10000000.0;
# boundsLow = -boundsLow*100000000.0


################################################################################################################
################################################################################################################
################################################################################################################


cgit  = 8;
maxit = 7;
alpha = 1e-1;
pcgTol = 1e-1;

HesPrec=getSSORCGRegularizationPreconditioner(1.0,1e-5,1000)

################################################# GIT VERSION OF JINV #################################################

regparams = [1.0,1.0,1.0,1e-6];
#### Use smoothness regularization 
# regfun(m, mref, M) = wdiffusionRegNodal(m, mref, M, Iact=Iact, C = regparams);
#### Use TV regularization 
regfun(m, mref, M) = wTVRegNodal(m, mref, M, Iact=Iact, C = regparams);
	
	
function dump(mc,Dc,iter,pInv,pMis)
	if resultsFilename!=""
		ntup = tuple((pInv.MInv.n+1)...);
		fullMc = slowSquaredToVelocity(reshape(Iact*pInv.modelfun(mc)[1] + sback,ntup))[1];
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

pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
#### Projected Gauss Newton
mc,Dc,flag = projGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);
#### Barrier Gauss Newton
# mc,Dc,flag = barrierGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump,epsilon = 0.1);

Dpred = Array(Array{Float64,2},length(pMisRFs))
for k = 1:length(pMisRFs)
	Dpred[k] = fetch(Dc[k]);
end
Dpred = arrangeRemoteCallDataIntoLocalData(Dpred);

if resultsFilename!=""
	Temp = splitext(resultsFilename);
	writedlm(string(Temp[1],"_predictedData",Temp[2]),convert(Array{Float16},Dpred));
	writedlm(string(Temp[1],"_recoveredModel",Temp[2]),convert(Array{Float16},reshape(Iact*modfun(mc)[1] + sback,tuple((pInv.MInv.n+1)...))));
end

#############################################################################################

