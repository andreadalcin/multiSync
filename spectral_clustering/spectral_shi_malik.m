function groups = spectral_shi_malik(CKSym,n,initial_datum)

MAXiter = 1000;
% Normalized spectral clustering according to Shi and Malik (2000)
% using Unnormalized Laplacian L = D - W
DN = diag(1./sqrt(sum(CKSym)+eps));
L = DN - CKSym;

[kerN,~] = eigs(L,DN,n,'smallestabs');

if initial_datum
    [~,ind] = sort(sum(CKSym),'descend');
    X0 = kerN(ind(1:n),:);
    groups = kmeans(kerN,n,'maxiter',MAXiter,'Start',X0,'EmptyAction','singleton','Distance','cityblock');
else
    REPlic = 100;
    groups = kmeans(kerN,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton','Distance','cityblock','Options',statset('UseParallel',1));
end

end