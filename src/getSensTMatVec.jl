export getSensTMatVec

function getSensTMatVec(v::Vector,m::Vector,pFor::EikonalInvParam)

   # extract pointers
    Q     = pFor.Sources
    P     = pFor.Receivers
    nsrc = size(Q,2); nrec = size(P,2)
	
	n = size(Q,1);

    # allocate space for result
    JTv  = zeros(n)
	JTvi = zeros(n)
	pEiks = pFor.eikonalParams;
	
	if length(v)!=nsrc*nrec
		warn("wrong size of vector v in EikonalInv::getSensTMatVec.jl");
	end
	
    for i=1:nsrc
		t = P*v[(i-1)*nrec + 1 : i*nrec];
		if pFor.useFilesForFields
			file = matopen(getFieldsFileName());
			pEiks[i].T1 		= read(file,string("T1_",i));
			pEiks[i].ordering = read(file,string("ordering_",i));
			pEiks[i].OP 		= read(file,string("OP_",i));
			close(file);
		end
		getSensTMatVecEik(t,JTvi,pEiks[i]);
		if pFor.useFilesForFields
			FactoredEikonalFastMarching.clear!(pEiks[i]);
		end
        JTv += JTvi;
    end
    return JTv;
end


