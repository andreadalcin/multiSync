function EdgeTable = build_PG_edgetable(Z,A,APG,clusters,R_patch,opts)

num_patches = length(clusters);
EdgeTable = empty_edgetable();

for i=1:num_patches
    for j=i+1:num_patches
        if APG(i,j) == 1
            M_i = R_patch{i};
            M_j = R_patch{j};
            
            A_filter = A;
            A_filter(setdiff(1:end,clusters{i}),:) = 0;
            A_filter(:,setdiff(1:end,clusters{j})) = 0;
            [row,col] = find(A_filter == 1);
            row = unique(row);
            col = unique(col);
            
            Mij_values = zeros(3,3,0);
            for r = 1:length(row)
                for c = 1:length(col)
                    idx_i = row(r);
                    idx_j = col(c);
                    if A(idx_i,idx_j) == 1
                        I1 = find(clusters{i} == idx_i);
                        I2 = find(clusters{j} == idx_j);
                        if I1 <= size(M_i,3) && I2 <= size(M_j,3)
                            X_i = M_i(:,:,I1);
                            X_j = M_j(:,:,I2);
                            if nnz(isnan(X_i(:))) == 0 && nnz(isnan(X_j(:))) == 0
                                Z_ij = Z(3*idx_i-2:3*idx_i,3*idx_j-2:3*idx_j);
                                Mij_values(:,:,end+1) = X_i \ Z_ij * X_j;
                            end
                        end
                    end
                end
            end
            
            if size(Mij_values,3) > opts.max_shared
                Mij_values = Mij_values(:,:,datasample(1:size(Mij_values,3),ceil(opts.max_shared/2)));
            end
            
            for k = 1:size(Mij_values,3)
                S = table(i,j,mat2cell(Mij_values(:,:,k),3,3),...
                    mat2cell(Mij_values(:,:,k)',3,3),...
                    'VariableNames',{'i','j','Z_ij','Z_ji'});
                EdgeTable = [EdgeTable; S];
            end
            
        end
    end
end

end