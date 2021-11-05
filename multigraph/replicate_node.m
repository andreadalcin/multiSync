function [S, corr] = replicate_node(E, i)

S = empty_edgetable();
corr = [];

% Iterate and look for edges
E_i = E(E.i == i,:);
last_idx = max([E.i; E.j]);

for j = i+1:max(E_i.j)
    E_ij = E_i(E_i.j == j,:);
    
    for idx = 1:size(E_ij,1)
        if idx == 1
            % Recycle the original vertex name in order not to leave blank
            % spaces
            k = i;
        else
            % i -> last_idx + idx
            k = last_idx + idx - 1;   
        end
        
        Z_ij = E_ij(idx,:).Z_ij;
        Z_ji = E_ij(idx,:).Z_ji;
        
        % Add new edge
        if k > j
            s = table(j,k,Z_ji,Z_ij,...
                'VariableNames',{'i','j','Z_ij','Z_ji'});
            S = [S; s];
        else
            s = table(k,j,Z_ij,Z_ji,...
                'VariableNames',{'i','j','Z_ij','Z_ji'});
            S = [S; s];
        end
        
        % Register a correspondence
        corr = unique([corr, k]);
    end
end
end