if nworkers()<=1
	addprocs(2);
end

using  jInv.Mesh
using  jInv.Utils
using  jInv.InverseSolve
using  EikonalInv
using  MAT

#############################################################################################################
modelDir = "../examples";


dataDir = pwd();
include("../drivers/prepareTravelTimeDataFiles.jl");
include("../drivers/plotModel.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");
include("../drivers/setupTravelTimeTomography.jl");

########################################################################################################
######################################## for SEG 256X512 ###############################################
#######################################################################################################
dim     = 2;
pad     = 0;
jump    = 40;
offset  = 256;
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[256,128],[0.0,13.5,0.0,4.2]);
dataFilenamePrefix = string(dataDir,"/DATA_SEG",tuple((Minv.n+1)...));


prepareTravelTimeDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,pad,jump,offset);
resultsFilenamePrefix = "";
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

cgit  = 2;
maxit = 1;
alpha = 1e+0;
pcgTol = 1e-1;

HesPrec=getSSORCGRegularizationPreconditioner(1.0,1e-5,1000)

################################################# GIT VERSION OF JINV #################################################

regparams = [1.0,1.0,1.0,1e-2];
regfun(m, mref, M) = wdiffusionRegNodal(m, mref, M, Iact=Iact, C = regparams);
	
pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
					 
mc,Dc,flag = projGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = []);

#############################################################################################
rm("DATA_SEG(256,128)_travelTime.dat");
rm("DATA_SEG(256,128)_rcvMap.dat");
rm("DATA_SEG(256,128)_srcMap.dat");
rm("DATA_SEG(256,128)_PARAM.mat");
rm("jInv.out");


