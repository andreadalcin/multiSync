function [dist_simple, dist_multi] = simple_multi_bench(N, ...
    num_multi, sigma)
% Runs a benchmark that compares the performance of a multigraph and of a
% simple graph obtained by combining rotations by averaging them.
% N is a 1x1 vector with the number of nodes to generate in the
% graph
% num_multi is a 1x1 vector containing the average multiplicity of 
% multi-edges
% sigma is a 1x1 vector with the standard deviation of the error to
% generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual implementation

% Generate a random graph
[~,Z,A] = random_graph(N);

% Solve the synchronization problem for the noiseless graph
synch_base = rotation_synch(Z,A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert simple graph
EdgeTable = graph2edgetable(Z,A);

% Add noise to the simple graph
E_1 = empty_edgetable();
E_2 = empty_edgetable();
p = 0.75;
n = size(A,1);
for i = 1:n-1
    for j = i+1:n
        if rand > (1-p)
            T = EdgeTable(EdgeTable.i == i & EdgeTable.j == j,:);
            
            R_tot = zeros(3);
            
            r = randi([1,2*num_multi], 1);
            for k = 1:r
                % Generate random rotation matrix
                R = rot_z(normrnd(0,sigma)) * rot_y(normrnd(0,sigma)) ...
                    * rot_x(normrnd(0,sigma));
                R_tot = R_tot + R;
                
                % Add multiedge with noise
                Z_ij = R * cell2mat(T.Z_ij);
                Z_ji = cell2mat(T.Z_ji) * R';
                S = table(i,j,mat2cell(Z_ij,3,3),mat2cell(Z_ji,3,3),...
                    'VariableNames',{'i','j','Z_ij','Z_ji'});
                E_1 = [E_1; S];
            end
            
            % Add a simple edge by averaging the rotation matrices
            % TO-DO: Improve rotation averaging. Govindu
            R_tot = R_tot ./ r;
            [U,~,V] = svd(R_tot);
            R_tot = U * diag([1,1,det(U*V')]) * V';
            
            Z_ij = R_tot * cell2mat(T.Z_ij);
            Z_ji = cell2mat(T.Z_ji) * R_tot';
            S = table(i,j,mat2cell(Z_ij,3,3),mat2cell(Z_ji,3,3),...
                'VariableNames',{'i','j','Z_ij','Z_ji'});
            E_2 = [E_2; S];
        end
    end
end

E_1 = sort_edges(E_1);
E_2 = sort_edges(E_2);


% Solve the synchronization problem for the multigraph
[E_1, R] = multigraph_expand(E_1);
[Z,A,C] = edgetable2graph(E_1,R);
synch_multi = constrained_synch(Z,A,C,R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the simple graph
[Z,A,~] = edgetable2graph(E_2,cell(1,N));
synch_simple = rotation_synch(Z,A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base = zeros(3,3,length(synch_base));
for k = 1:length(synch_base)
    base(:,:,k) = synch_base{k};
end

simple = zeros(3,3,length(synch_simple));
for k = 1:length(synch_simple)
    simple(:,:,k) = synch_simple{k};
end

multi = zeros(3,3,length(synch_multi));
for k = 1:length(synch_multi)
    multi(:,:,k) = synch_multi{k};
end

dist_simple = error_R(simple,base);
dist_multi = error_R(multi,base);

end