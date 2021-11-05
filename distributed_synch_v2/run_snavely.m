function results = run_snavely(Z,A,R_gt,cameras_gt,i_cc,i_gt,opts)

% Synchronization on patches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create patches from graph
n = size(A,1);
clusters = create_graph_patches(A,opts);

% SO(3) synchronization for all patches
num_patches = length(clusters);
R_patch = cell(num_patches,1);

A_clustered = cell(num_patches);
Z_clustered = cell(num_patches);
irls = cell(num_patches);
opts.max_size_coefficient = 10;
for k = 1:num_patches
    A_clustered{k} = A(clusters{k},:);
    Z_clustered{k} = Z;
    irls{k} = opts.irls;
end

tic
parfor k = 1:num_patches
    % Get subset of adjacency matrix corresponding to k-th cluster
    A_patch = A_clustered{k};
    A_patch = A_patch(:,clusters{k});
    
    % Get subset of Z corresponding to k-th cluster
    index_patch = sort([3*clusters{k}-2; 3*clusters{k}-1; 3*clusters{k}],...
        'ascend');
    Z = Z_clustered{k};
    Z_patch = Z(index_patch,:);
    Z_patch = Z_patch(:,index_patch);
    
    % Solve synchronization for k-th cluster
    % IRLS
    % R_patch{k} = SO3_eig_IRLS(Z_patch,A_patch,30,1e-5,'cauchy',1,0,0);
    % EIGS
    % R_patch{k} = SO3_EIG(Z_patch,A_patch);
    % MPLS
    R_patch{k} = rotation_synch_MPLS(Z_patch,A_patch,opts);
    % L1-IRLS
%     fprintf("Testing: L1-IRLS (Govindu)...\n");
%     [I,J] = find(A_patch);
%     n_pairs = length(I);
%     R_pairwise = zeros(3,3,n_pairs);
%     for p = 1:n_pairs
%         j = I(p);
%         i = J(p);
%         R_pairwise(:,:,p) = Z_patch(3*i-2:3*i,3*j-2:3*j);
%     end
%     R_patch{k} = AverageSO3Graph(R_pairwise,[I J]');
    % R-GoDec
%         power = 1;
%         iter_max = 100;
%         rank = 3;
%         n_rotations = nnz(A_patch);
%         lambda = 0.02* sqrt(2*log(9*n_rotations));
%         ncams = size(A_patch,1);
%         [~,ind]=max(sum(A_patch));
%         [G_godec,~,~] = R_GoDec_mixed(Z_patch,rank,lambda,power,iter_max,1e-5);
%         R_patch{k} = RijtoRi(G_godec,ncams,ind);
end
t1 = toc;
fprintf("Time for rotation synch on clusters: %.2f\n", t1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build patch graph and solve synchronization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
A_PG = build_PG_A(A, clusters);
X_PG = build_PG_X(Z, A, A_PG, clusters, R_patch, opts);
edgetable_PG = build_PG_edgetable(Z, A, A_PG, clusters, R_patch, opts);
t2 = toc;
fprintf("Time patch graph build: %.2f\n", t2);

%==========================================================================
% Solve ambiguity for edge aggregation

% Synchronization on patch graph
tic
R_PG_aggr = rotation_synch_MPLS(X_PG,A_PG,opts);
t3 = toc;

% Apply rotation from patch graph to rotations from patches
R_patch_aggr = cell(num_patches,1);
for k = 1:num_patches
    for i = 1:length(clusters{k})
        R_patch_aggr{k}(:,:,i) = R_patch{k}(:,:,i) * (R_PG_aggr(:,:,k));
    end
end

% Merge rotations
R_aggr = zeros(3,3,n);
for i = 1:n
    R_values = [];
    count = 0;
    
    for k = 1:num_patches
        ind = find(clusters{k} == i,1);
        if ~isempty(ind)
            count = count+1;
            R_values(:,:,count) = R_patch_aggr{k}(:,:,ind);
            break;
        end
    end
    
    % Perform single averaging
    R_aggr(:,:,i) = R_values;
end


fprintf("Time - Edge averaging: %.2f\n", t1 + t2 + t3);


%==========================================================================
% Solve ambiguity for multigraph synch

% % Synchronization on patch graph w/ multigraph
[E_multi, R_multi, t] = multigraph_expand(edgetable_PG);
tic
[Z_multi,A_multi,C_multi] = edgetable2graph(E_multi,R_multi);
R_PG_multi = constrained_synch_IRLS(Z_multi,A_multi,C_multi,R_multi,30,1e-5,'cauchy',1);
t4 = toc;

% Apply rotation from patch graph to rotations from patches
R_patch_multi = cell(num_patches,1);
for k = 1:num_patches
    for i = 1:length(clusters{k})
        R_patch_multi{k}(:,:,i) = R_patch{k}(:,:,i) * R_PG_multi{k};
    end
end

% Merge rotations
R_multi = zeros(3,3,n);
for i = 1:n
    R_values = [];
    count = 0;

    for k = 1:num_patches
        ind = find(clusters{k} == i,1);
        if ~isempty(ind)
            count = count+1;
            R_values(:,:,count) = R_patch_multi{k}(:,:,ind);
            break;
        end
    end

    % Perform single averaging
    R_multi(:,:,i) = R_values;
end


fprintf("Time - Multigraph synch: %.2f\n", t1 + t2 + t + t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_standard = error_R(R_aggr(:,:,i_cc),R_gt(:,:,cameras_gt(i_gt)+1));
err_multi = error_R(R_multi(:,:,i_cc),R_gt(:,:,cameras_gt(i_gt)+1));

results.err_standard = err_standard;
results.err_multi = err_multi;

end