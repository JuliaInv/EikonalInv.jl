using jInv.Mesh
using EikonalInv
using jInv.ForwardShare
using jInv.Utils

using Base.Test



function testSensMatVec(pFor,m0)

	println("\t\tchecking derivatives")

	dobs, = getData(vec(m0),pFor)


	function testFun(m,v=[])
		dobs, = getData(vec(m),pFor)
		if isempty(v)
			return vec(dobs)
		else
			JTv = getSensMatVec(vec(v),vec(m),pFor)
			return vec(dobs), JTv
		end
	end

	chkDer, = jInv.Utils.checkDerivative(testFun,vec(m0))
	@test chkDer

	(nd,nm) = getSensMatSize(pFor)
	v0 = randn(nm)
	w0 = randn(nd)

	J  = getSensMat(vec(m0),pFor)
	J1 = J*v0;
	J2 = J'*w0;

	J1mf = getSensMatVec(v0,vec(m0),pFor)
	J2mf = getSensTMatVec(w0,vec(m0),pFor)

	println("\t\ttest matrix free sensitivities")
	@test norm(J1mf-J1)/norm(J1) < 1e-10
	@test norm(J2mf-J2)/norm(J2) < 1e-10

	println("\t\tadjoint test")
	t1 = dot(w0,J1mf)
	t2 = dot(v0,J2mf)


	p = abs(t1-t2)/abs(t1);
	@test p < 1e-10
	if p > 1e-10
		println("Inconsistency in sensitivity: ",p);
	end

end
############################################################################################
# get 2D mesh with uneven number of cells
domain = [1. 4.0 .4 4.2]
n      = [17;16;]
M      = getRegularMesh(domain,n)
# sources / receivers on top edge
idx    = reshape(collect(1:prod(n+1)),tuple(n+1...))
ib     = idx[1:2:end,1];
I      = speye(prod(n+1))
Q      = I[:,vec(ib)]
R      = copy(Q)


println("2D - test low order scheme")
pFor = getEikonalInvParam(M,Q,R,false)
m0   = rand(prod(n+1))+10.0

try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for low-order scheme and mesh size n=$n FAILED!")
end

############################################################################################
println("2D - test high order scheme")
pFor = getEikonalInvParam(M,Q,R,true)
m0   = rand(prod(n+1))+10.0
try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for higher-order scheme and mesh size n=$n FAILED!")
end

n      = [16;32;]
M      = getRegularMesh(domain,n)
# sources / receivers on top edge
idx    = reshape(collect(1:prod(n+1)),tuple(n+1...))
ib     = idx[1:2:end,1];
I      = speye(prod(n+1))
Q      = I[:,vec(ib)]
R      = copy(Q)
println("2D - test low order scheme")
pFor = getEikonalInvParam(M,Q,R,false)
m0   = rand(prod(n+1))+10.0
try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for lower-order scheme and mesh size n=$n FAILED!")
end

println("2D - test high order scheme")
pFor = getEikonalInvParam(M,Q,R,true)
m0   = rand(prod(n+1))+10.0
try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for higher-order scheme and mesh size n=$n FAILED!")
end
############################################################################################

# get 3D mesh
domain = [1. 4.0 .4 4.2 0 1.2]
n      = [16;16;8]
M      = getRegularMesh(domain,n)
# sources / receivers on top edge
idx    = reshape(collect(1:prod(n+1)),tuple(n+1...))
ib     = idx[1:6:end,1,1:4:end];
I      = speye(prod(n+1))
Q      = I[:,vec(ib)]
R      = copy(Q)
println("3D - test low order scheme")
pFor = getEikonalInvParam(M,Q,R,false)
m0   = rand(prod(n+1))+10.0
try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for low-order scheme and mesh size n=$n FAILED!")
end

println("3D - test high order scheme")
pFor = getEikonalInvParam(M,Q,R,true)
m0   = rand(prod(n+1))+10.0
try
	testSensMatVec(pFor,vec(m0))
catch
	warn("Sensitivity for higher-order scheme and mesh size n=$n FAILED!")
end
