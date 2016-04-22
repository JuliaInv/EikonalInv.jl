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
include("../drivers/invertTravelTimeTomography.jl");

NumIter = 1;

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
(mc,Dpred) = invertTravelTimeTomographyBoundModel(m,dataFilenamePrefix, resultsFilenamePrefix,NumIter);
#############################################################################################
rm("DATA_SEG(256,128)_travelTime.dat");
rm("DATA_SEG(256,128)_rcvMap.dat");
rm("DATA_SEG(256,128)_srcMap.dat");
rm("DATA_SEG(256,128)_PARAM.mat");


