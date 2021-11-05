function APG = build_PG_A(A,clusters)

num_patches = length(clusters);
APG = eye(num_patches);

for i = 1:num_patches
    for j = i+1:num_patches
        A_ij = A(clusters{i},:);
        A_ij = A_ij(:,clusters{j});
        if nnz(A_ij)
            APG(i,j) = 1;
            APG(j,i) = 1;
        end
    end
end

APG = sparse(APG);

end