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
		getSensTMatVecEik(t,JTvi,pEiks[i]);
        JTv += JTvi;
    end
    return JTv;
end


