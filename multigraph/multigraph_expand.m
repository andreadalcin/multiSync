function [E, R, t] = multigraph_expand(EdgeTable)

n = max([EdgeTable.i; EdgeTable.j]);
E = sort_edges(EdgeTable);
R = cell(1,n);  % Store correspondences between nodes and replicas
Q = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimized node expansion
tic
for i = 1:n-1
    count = count_multiedges(E,i,n);
    if count > 1
        % Replicate node
        [S, corr] = replicate_node(E, i);
        
        % Delete starting node before adding in order not to overwrite
        E(E.i == i,:) = [];
        
        % Add new edges to the graph and save indexes of node replicas
        E = [E; S];
        R{i} = corr;
    elseif count == 1
        Q = [Q; i];
    end
end

for i = 1:size(Q,1)
    if count_multiedges(E,Q(i),n) >= 1
        % Replicate node
        [S, corr] = replicate_node(E, Q(i));
        
        % Delete starting node before adding in order not to overwrite
        E(E.i == Q(i),:) = [];
        
        % Add new edges to the graph and save indexes of node replicas
        E = [E; S];
        R{Q(i)} = corr;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add identity constraints
for k = 1:n
    if ~isempty(R{k})
        pairs = nchoosek(R{k}, 2);
        for h = 1:size(pairs,1)
            S = table(pairs(h,1), pairs(h,2), mat2cell(eye(3),3,3), ...
                mat2cell(eye(3),3,3),'VariableNames',{'i','j','Z_ij','Z_ji'});
            E = [E; S];
        end
    end
end

E = sort_edges(E);
t = toc;

end