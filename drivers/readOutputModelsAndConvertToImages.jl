using PyPlot

A = readdir();

for k=1:length(A)
	fname = A[k];
	if fname[end-3:end] == ".dat"
		B = readdlm(fname);
		figure(999)
		imshow(B');
		colorbar();
		title(fname[1:end-4],fontsize = 11);
		savefig(string(fname[1:end-4],".png"));
		close(999)
	end
end