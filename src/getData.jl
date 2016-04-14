export getData

function getData(m,pFor::EikonalInvParam,doClear::Bool=false)

    # extract pointers
    Mesh   	= pFor.Mesh
    Q     	= pFor.Sources
    P     	= pFor.Receivers
	n_nodes = pFor.Mesh.n+1;
    nrec  	= size(P,2) 
    nsrc  	= size(Q,2)
    
    # allocate space for data and fields
    D  = zeros(nrec,nsrc)
    pEik = Array(EikonalParam,nsrc);
    
	pMem = getEikonalTempMemory(n_nodes);
	ntup = tuple(n_nodes...);
	T = zeros(Float64,ntup) # V is a temporary array of Float64. 
    for k=1:nsrc
		src_k_loc = Q.rowval[k];
		if Mesh.dim==2
			src = zeros(Int64,2);
			cs2loc(src,src_k_loc,n_nodes);
		else
			src = zeros(Int64,3);
			cs2loc3D(src,src_k_loc,n_nodes);
		end
		m = reshape(m,ntup);
		pEik[k] = getEikonalParam(Mesh,m,src,pFor.HO);
		pEik[k].T1 = T; # Here, we set T (of Float64) for the calculation, but reuse the memory.
		solveFastMarchingUpwindGrad(pEik[k],pMem);
		pEik[k].T1 = zeros(Float32,ntup); # For the sensitivities, we don't really need a high precision T.
		pEik[k].T1[:] = T;
		if Mesh.dim==2
			selfMultiplyWithAnalyticSolution2D(n_nodes,Mesh.h,src,T);
		else
			selfMultiplyWithAnalyticSolution3D(n_nodes,Mesh.h,src,T)
		end
		D[:,k] = P'*T[:];
		if doClear
			FactoredEikonalFastMarching.clear!(pEik[k]);
		end
	end
	pFor.eikonalParams = pEik;
    return D,pFor
end