export getBoundModel,getBoundModelInv
function getBoundModel(m,a,b) 
	d = (b-a)./2.0;
	mid = (b+a)./2.0;
	v = (1.0./d).*(m - mid);
	s = d.*(tanh(v)+1) + a;
	dsdm = spdiagm(sech(v).^2)
	return s, dsdm
end

function getBoundModelInv(s,a,b)
	d = (b-a)./2.0;
	mid = (b+a)./2.0;
    m = atanh((s-a).*(1.0./d)-1) .* d + mid;
	return m
end