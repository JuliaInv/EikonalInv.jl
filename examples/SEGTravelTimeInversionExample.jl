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


modelDir = pwd();
include("../drivers/prepareTravelTimeDataFiles.jl");
include("../drivers/plotModel.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");
include("../drivers/invertTravelTimeTomography.jl");

NumIter = 10;

########################################################################################################
######################################## for SEG 256X512 ###############################################
#######################################################################################################
dim     = 2;
pad     = 0;
jump    = 20;
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
resultsFilenamePrefix = "";
# this filename can be used to save the model during the iterations
# resultsFilenamePrefix = string(resultsDir,"/travelTimeInvSEG"); 
(mc,Dpred) = invertTravelTimeTomographyBoundModel(m,dataFilenamePrefix, resultsFilenamePrefix,NumIter);

if plotting
	figure()
	plotModel(mc,true,true,Minv,pad,limits);
end

#############################################################################################

