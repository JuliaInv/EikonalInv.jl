export expandModelNearest, getSimilarLinearModel, addAbsorbingLayer,cutAbsorbingLayer

function expandModelNearest(m,n,ntarget)
if length(size(m))==2
	mnew = zeros(Float64,ntarget[1],ntarget[2]);
	for j=1:ntarget[2]
		for i=1:ntarget[1]
			jorig = convert(Int64,ceil(j/ntarget[2]*n[2]));
			iorig = convert(Int64,ceil(i/ntarget[1]*n[1]));
			mnew[i,j] = m[iorig,jorig];
		end
	end
elseif length(size(m))==3
	mnew = zeros(Float64,ntarget[1],ntarget[2],ntarget[3]);
	for k=1:ntarget[3]
		for j=1:ntarget[2]
			for i=1:ntarget[1]
				korig = convert(Int64,ceil(k/ntarget[3]*n[3]));
				jorig = convert(Int64,ceil(j/ntarget[2]*n[2]));
				iorig = convert(Int64,ceil(i/ntarget[1]*n[1]));
				mnew[i,j,k] = m[iorig,jorig,korig];
			end
		end
	end
end
return mnew
end

function getSimilarLinearModel(m::Array{Float64})
if length(size(m))==2
	(nx,nz) = size(m);
	m_vel = 1./sqrt(m);
	mtop = m_vel[:,5:20];
	mtop = mean(mtop[:]);
	mbottom = m_vel[:,end-30:end-10];
	mbottom = mean(mbottom[:]);
	m_vel = ones(nx)*linspace(mtop,mbottom,nz)';
	mref = 1./(m_vel.^2);
	
	
	# println("~~~~~~~~~~ getSimilarLinearModel ~~~~~~~~~~~~~~~~~~~~~~~~")
	# println(mtop);println(mbottom);
	# X,Y = ndgrid(0:size(m,1)-1,0:size(m,2)-1);
	# y1 = 1./sqrt(m[1,1]);
	# y2 = 1./sqrt(m[1,end]);
	# a = (y2-y1)/size(m,2);
	# vel0 = y1+a*Y;
	# mref = ((1./vel0).^2);
	
	
	
	
	# vel_left = reshape(1./sqrt(m[1,:]),size(m,2));
	# vel_right = reshape(1./sqrt(m[end,:]),size(m,2));
	# z = 0:size(m,2)-1;
	# a = (vel_left[1] + vel_right[1]) / 2.0;
	# w = 1-(z ./ maximum(z));
	# w = (z/maximum(z));
	# b = (dot(vel_left - a,w.*z) + dot(vel_right - a,w.*z) )/ (2*dot(z,w.*z))
	# profile = 1./((b*z + a).^2);
	# mref = ones(size(m,1))*profile'
	
	
	# vel_left = reshape(1./sqrt(m[1,:]),size(m,2));
	# vel_right = reshape(1./sqrt(m[end,:]),size(m,2));
	# z = 0:size(m,2)-1;
	# vel_top = (vel_left[1] + vel_right[1]) / 2.0;
	# vel_bot = 1./sqrt(m[div(size(m,1),2),end]);
	# b = (vel_bot - vel_top) / size(m,2)
	# profile = 1./((b*z + vel_top).^2);
	# mref = ones(size(m,1))*profile'
	
	
	
	
	# println("~~~~~~~~~~ getSimilarLinearModel Old ~~~~~~~~~~~~~~~~~~~~~~")
	# println(y1);println(y2);
	
   

elseif length(size(m))==3
	(nx,ny,nz) = size(m);
	m_vel = 1./sqrt(m);
	mtop = m_vel[:,:,5:15];
	mtop = mean(mtop[:]);
	mbottom = m_vel[:,:,end-30:end-10];
	mbottom = mean(mbottom[:]);
	lin = linspace(mtop,mbottom,nz);
	m_vel = copy(m);
	Oplane = ones(nx,ny);
	for k=1:nz
		m_vel[:,:,k] = lin[k]*Oplane;
	end
	mref = 1./(m_vel.^2);
else
	error("Unhandled Dimensions");
end
return mref;
end


function addAbsorbingLayer2D(m::Array{Float64},pad::Int64)
if pad<=0
	return m;
end
mnew = zeros(size(m,1)+2*pad,size(m,2)+pad);
mnew[pad+1:end-pad,1:end-pad] = m;
mnew[1:pad,1:end-pad] = repmat(m[1,:],pad,1);
mnew[end-pad+1:end,1:end-pad] = repmat(m[end,:],pad,1);
mnew[:,end-pad+1:end] = repmat(mnew[:,end-pad],1,pad);

return mnew;
end


function addAbsorbingLayer(m::Array{Float64},Msh::RegularMesh,pad::Int64)
if pad<=0
	return m,Msh;
end
Omega = Msh.domain;

if length(size(m))==2
	mnew = addAbsorbingLayer2D(m,pad);
	MshNew = getRegularMesh([Omega[1],Omega[2] + 2*pad*Msh.h[1],Omega[3],Omega[4]+pad*Msh.h[2]],collect(size(mnew))-1);
elseif length(size(m))==3
	mnew = zeros(size(m,1)+2*pad,size(m,2)+2*pad,size(m,3)+pad);
	mnew[pad+1:end-pad,pad+1:end-pad,1:end-pad] = m;
	
	extendedPlane1 = addAbsorbingLayer2D(reshape(m[1,:,:],size(m,2),size(m,3)),pad);
	extendedPlaneEnd = addAbsorbingLayer2D(reshape(m[end,:,:],size(m,2),size(m,3)),pad);
	
	for k=1:pad
		mnew[k,:,:] = extendedPlane1;
		mnew[end-k+1,:,:] = extendedPlaneEnd;
		mnew[pad+1:end-pad,end-k+1,1:end-pad] = m[:,end,:];
		mnew[pad+1:end-pad,k,1:end-pad] = m[:,1,:];
	end
	t = mnew[:,:,end-pad];
	for k=1:pad
		mnew[:,:,end-pad+k] = t;
	end
	MshNew = getRegularMesh([Omega[1],Omega[2] + 2*pad*Msh.h[1],Omega[3],Omega[4] + 2*pad*Msh.h[2],Omega[5],Omega[6]+pad*Msh.h[2]],collect(size(mnew))-1);
end

return mnew,MshNew;
end


function cutAbsorbingLayer(m::Array{Float64},Msh,pad::Int64)
if pad<=0
	return m,Msh;
end
Omega = Msh.domain;
if length(size(m))==2
	mnew = m[pad+1:end-pad,1:end-pad];
	nnew = collect(size(mnew))-1;
	OmegaNew = [Omega[1],Omega[2] - 2*pad*Msh.h[1],Omega[3],Omega[4]-pad*Msh.h[2]];
	MshNew = getRegularMesh(OmegaNew,nnew);
elseif length(size(m))==3
	error("Not implemented")
end

return mnew,MshNew;
end