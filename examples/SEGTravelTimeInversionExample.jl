
using jInv.Mesh
using jInv.Utils
using jInv.InverseSolve

using EikonalInv
using MAT

#############################################################################################################
plotting = false;
if plotting
	using jInvVis
	using  PyPlot
	close("all");
end
#############################################################################################################
dataDir = pwd();
resultsDir = pwd();

###############################################################################################################
# this filename can be used to save the model during the iterations. If set to "", no images or files will be saved.
# resultsFilenamePrefix = "";
resultsFilenamePrefix = string(resultsDir,"/shortTravelTimeInvSEG"); 
################################################################################################################


modelDir = pwd();
include("../drivers/prepareTravelTimeDataFiles.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");
include("../drivers/setupTravelTimeTomography.jl");

########################################################################################################
######################################## for SEG 256X512 ###############################################
#######################################################################################################
dim     = 2;
pad     = 0;
jump    = 5;
newSize = [256,128];
offset  = ceil(Int64,(newSize[1]*(6/13.5)));;
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[0.0,13.5,0.0,4.2],newSize,1.752,2.9);
dataFilenamePrefix = string(dataDir,"/shortDATA_SEG",tuple((Minv.n+1)...));

limits = (minimum(m),maximum(m));
	
if plotting
	figure()
	plotModel(m,true,Minv,pad,limits);
	figure()
	plotModel(mref,true,Minv,pad,limits);
end

prepareTravelTimeDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,pad,jump,offset);

(Q,P,pMisRFs,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename,Dobs) = setupTravelTimeTomography(m,dataFilenamePrefix, resultsFilenamePrefix);

########################################################################################################
##### Set up Inversion #################################################################################
########################################################################################################
############### Inversion for the velocity:
# modfun = velocityToSlowSquared;
################################################################################################################
############### USE A BOUND MODEL For velocity #################################################################
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
# boundsLow = -boundsLow*100000000.0;
################################################################################################################
################ Inversion for the squared slowness: ###########################################################
################################################################################################################

mref 		= velocityToSlowSquared(mref)[1];
t    		= copy(boundsLow*0.9);
boundsLow 	= velocityToSlowSquared(boundsHigh*1.1)[1];
boundsHigh 	= velocityToSlowSquared(t)[1]; t = 0;
###################################################
a = minimum(boundsLow);
b = maximum(boundsHigh);
maxStep=0.1*maximum(boundsHigh);


############### Standard bounded optimization:
modfun 		= identityMod;

############### USE A BOUND MODEL For slow squared #############################################################

# modfun(m) = getBoundModel(m,a,b);
# mref = getBoundModelInv(mref,a,b);
# boundsHigh = boundsHigh*10000000.0;
# boundsLow = -boundsLow*100000000.0


################################################################################################################
################################################################################################################
################################################################################################################


cgit  = 8;
maxit = 10;
alpha = 5e-1;
pcgTol = 1e-1;

HesPrec=getSSORCGRegularizationPreconditioner(1.0,1e-5,100)

################################################# GIT VERSION OF JINV #################################################

regparams = [1.0,1.0,1.0,1e-6];
#### Use smoothness regularization 
regfun(m, mref, M) = wdiffusionRegNodal(m, mref, M, Iact=Iact, C = regparams);
#### Use TV regularization 
# regfun(m, mref, M) = wTVRegNodal(m, mref, M, Iact=Iact, C = regparams);
	
	
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
			# plotModel(fullMc,false,[],0,limits,splitdir(Temp)[2]); # no information on axes
			plotModel(fullMc,true,Minv,pad,limits,splitdir(Temp)[2]); # with domain information on axes
		end
	end
end	

pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);

D0RFs, = computeMisfit(pInv.modelfun(mref[:])[1],pMisRFs,false);
					 
#### Projected Gauss Newton
mc,Dc,flag,his = projGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);
#### Barrier Gauss Newton
# mc,Dc,flag,his = barrierGNCG(copy(mref[:]),pInv,pMisRFs,indCredit = [],dumpResults = dump);


if plotting
	figure()
	Fvals = his.F;
	Jvals = his.Jc;
	numiter = length(Fvals);
	semilogy(0:numiter-1,Fvals,"-*");
	xlabel("Gauss Newton iterations")
	ylabel("Misfit value")
end

D0 = Array(Array{Float64,2},length(pMisRFs))
for k = 1:length(pMisRFs)
	D0[k] = fetch(D0RFs[k]);
end
D0 = arrangeRemoteCallDataIntoLocalData(D0); D0RFs = 0;



Dpred = Array(Array{Float64,2},length(pMisRFs))
for k = 1:length(pMisRFs)
	Dpred[k] = fetch(Dc[k]);
end
Dpred = arrangeRemoteCallDataIntoLocalData(Dpred); Dc = 0;

if plotting
	figure()
	imshow(Dpred');
	colorbar();
	xlabel("Sources")
	ylabel("Receivers")
	figure()
	imshow(D0');
	colorbar();
	xlabel("Sources")
	ylabel("Receivers")
	figure()
	imshow(abs(Dpred - Dobs)');
	colorbar();
	xlabel("Sources")
	ylabel("Receivers")
	figure()
	imshow(abs(D0 - Dobs)');
	colorbar();
	xlabel("Sources")
	ylabel("Receivers")
end

if resultsFilename!=""
	Temp = splitext(resultsFilename);
	writedlm(string(Temp[1],"_predictedData",Temp[2]),convert(Array{Float16},Dpred));
	writedlm(string(Temp[1],"_recoveredModel",Temp[2]),convert(Array{Float16},reshape(Iact*modfun(mc)[1] + sback,tuple((pInv.MInv.n+1)...))));
end

#############################################################################################

