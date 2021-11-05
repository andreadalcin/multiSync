function groups = spectral_cluster(CKSym,n,initial_datum)
N = size(CKSym,1);
MAXiter = 1000;
DN = diag( 1./sqrt(sum(CKSym)+eps) );
LapN = speye(N) - DN * CKSym * DN;
[~,~,kerN] = svds(LapN,n,0);
for i = 1:N
    kerN(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
end

if initial_datum
    [~,ind]=sort(sum(CKSym),'descend');
    X0=kerN(ind(1:n),:);
    groups = kmeans(kerN,n,'maxiter',MAXiter,'Start',X0,'EmptyAction','singleton','Distance','cityblock');
else
    REPlic = 100;
    groups = kmeans(kerN,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton','Distance','cityblock','Options',statset('UseParallel',1));
end

end