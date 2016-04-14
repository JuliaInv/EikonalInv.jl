export getSensMatVec

function getSensMatVec(v::Vector,m::Vector,pFor::EikonalInvParam)

    # extract pointers
    Q     = pFor.Sources
    P     = pFor.Receivers

    nsrc = size(Q,2)
    
    # allocate space for matvec product
    Jv   = zeros(size(P,2),nsrc)
	
	if length(v)!=length(m)
		warn("wrong size of vector v in EikonalInv::getSensMatVec.jl");
	end
	
    
    # derivative of mass matrix
	Tlin = zeros(size(v));
	
	if nsrc != length(pFor.eikonalParams)
		error("EikonalInv: Calculating derivatives without having fields... Run getData() to pFor");
	end
	
	pEiks = pFor.eikonalParams;
    for i=1:nsrc
		pEik_i = pEiks[i];
		getSensMatVecEik(v,Tlin,pEik_i);
		Jv[:,i] = P'*Tlin;
    end
    return vec(Jv)
end