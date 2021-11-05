function  B=adj2inc(A)
nV  = size(A,1);
[vNodes1,vNodes2] = find(tril(A,-1));
nE = length(vNodes1);
vOnes = ones(nE,1);
vEdgesidx = 1:nE;

B= sparse([vNodes1; vNodes2],...
    [vEdgesidx, vEdgesidx]',...
    [vOnes; -vOnes],...
    nV,nE);


end

