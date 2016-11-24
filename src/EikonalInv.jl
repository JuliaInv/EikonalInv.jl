module EikonalInv

using jInv.Mesh
using jInv.Utils
using FactoredEikonalFastMarching
using MAT


import jInv.ForwardShare.getData
import jInv.ForwardShare.getSensTMatVec
import jInv.ForwardShare.getSensMatVec

import jInv.ForwardShare.ForwardProbType

export EikonalInvParam
export getEikonalInvParam

function getFieldsFileName()
	tfilename = string("tempEikFields_worker",myid(),".mat");
end

"""
type EikonalInvParam

A type to be used with jInv for solving the inverse eikonal equation (e.g., travel time tomography).

Fields:

	Mesh 			:: RegularMesh
	Sources 		:: SparseMatrixCSC  		  - The matrix of point source locations of size N x num_sources
	Receivers 		:: SparseMatrixCSC  		- The matrix of point reciever locations of size N x num_recievers
	HO          	::Bool              		- flag for higher-order discretization
	eikonalParams	::Array{EikonalParam}   - An array of forward problem parameters. See FactoredEikonalFastMarching for EikonalParam.
  useFilesForFields :: Bool             - flag for writing fields to disk

Constructor:

	getEikonalInvParam
"""
type EikonalInvParam <: ForwardProbType
	Mesh      		  :: RegularMesh
  Sources			  :: SparseMatrixCSC
  Receivers		  :: SparseMatrixCSC
	HO        		  :: Bool
	eikonalParams	  :: Array{EikonalParam}
	useFilesForFields :: Bool
end

import jInv.ForwardShare.getSensMatSize
getSensMatSize(pFor::EikonalInvParam) = (size(pFor.Receivers,2)*size(pFor.Sources,2),prod(pFor.Mesh.n+1))

"""
function EikonalInv.getEikonalInvParam

Constructor for EikonalInvParam

Input:

	Mesh::RegularMesh
	Sources::SparseMatrixCSC      - The matrix of point source locations of size N x num_sources
	Receivers::SparseMatrixCSC    - The matrix of point reciever locations of size N x num_recievers
	HO::Bool=false							  - flag for higher-order discretization
	useFilesForFields::Bool=false - flag for writing fields to disk
"""
function getEikonalInvParam(Mesh::RegularMesh,Sources::SparseMatrixCSC,Receivers::SparseMatrixCSC,HO::Bool=false,useFilesForFields::Bool = false)
	## This function does not use the parallel mechanism of jInv.
	return EikonalInvParam(Mesh,Sources,Receivers,HO,Array(EikonalParam,0),useFilesForFields);
end


function getEikonalInvParam(Mesh::RegularMesh,Sources::SparseMatrixCSC,Receivers::SparseMatrixCSC,HO::Bool,numWorkers::Int64,useFilesForFields::Bool = false)
	## This function does use the parallel mechanism of jInv (i.e., returns a RemoteChannel), even if numWorkers=1.
	if numWorkers > nworkers()
		numWorkers = nworkers();
	end
	SourcesSubInd = Array(Array{Int64,1},numWorkers);
	ActualWorkers = workers();
	if numWorkers < nworkers()
		ActualWorkers = ActualWorkers[1:numWorkers];
	end
	continuationDivision = [1;numWorkers+1];

	pFor   = Array{RemoteChannel}(numWorkers)
	i = 1; nextidx() = (idx=i; i+=1; idx)

	nsrc  = size(Sources,2);
	# send out jobs
	@sync begin
		for p=workers()
			@async begin
				while true
					idx = nextidx()
					if idx > numWorkers
						break
					end
					I_p = getSourcesIndicesOfKthWorker(numWorkers,idx,nsrc);
					SourcesSubInd[idx] = I_p;
					# find src and rec on mesh
					pFor[idx]  = initRemoteChannel(getEikonalInvParam,p, Mesh, Sources[:,I_p], Receivers, HO,useFilesForFields);
					wait(pFor[idx])
				end
			end
		end
	end
	return pFor,continuationDivision,SourcesSubInd# Array of Remote Refs
end

import jInv.Utils.clear!
function clear!(pFor::EikonalInvParam)
for k = 1:length(pFor.eikonalParams)
	clear!(pFor.eikonalParams[k]);
end
return pFor;
end

include("./SeismicUtils/SeismicUtils.jl")
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")


end
