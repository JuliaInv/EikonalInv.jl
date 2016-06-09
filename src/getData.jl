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
    
	
	if pFor.useFilesForFields
		tfilename = getFieldsFileName();
		tfile = matopen(tfilename, "w");
	end
	
	pMem = getEikonalTempMemory(n_nodes);
	ntup = tuple(n_nodes...);
	T = zeros(Float64,ntup) # V is a temporary array of Float64. 
    T1_temp = zeros(Float32,ntup);
	ordering_temp = zeros(FactoredEikonalFastMarching.EikIdxType,prod(n_nodes));
	OP_temp 	  = zeros(Int8,ntup)
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
		pEik[k].ordering = ordering_temp;
		pEik[k].OP		 = OP_temp;
		solveFastMarchingUpwindGrad(pEik[k],pMem);
		T1_temp[:] = T;
		if Mesh.dim==2
			selfMultiplyWithAnalyticSolution2D(n_nodes,Mesh.h,src,T);
		else
			selfMultiplyWithAnalyticSolution3D(n_nodes,Mesh.h,src,T)
		end
		D[:,k] = P'*T[:];
		if doClear
			FactoredEikonalFastMarching.clear!(pEik[k]);
		elseif pFor.useFilesForFields
			write(tfile,string("T1_",k),pEik[k].T1);
			write(tfile,string("ordering_",k),pEik[k].ordering);
			write(tfile,string("OP_",k),pEik[k].OP);
			FactoredEikonalFastMarching.clear!(pEik[k]);
		else
			# Only here we actually save the data.
			pEik[k].T1 = copy(T1_temp);
			pEik[k].ordering = copy(ordering_temp);
			pEik[k].OP = copy(OP_temp);
		end
	end
	if pFor.useFilesForFields
		close(tfile);
		gc()
	end
	pFor.eikonalParams = pEik;
    return D,pFor
end