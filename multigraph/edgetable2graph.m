% Converts an edge table to a graph
% E is the edge table
% R is an array of cell with node replicas correspondences
function [Z,A,C] = edgetable2graph(E,R)
n = max([E.i; E.j]);

% Check that the edge table has no multiedges
for i=min([E.i;E.j]):n-1
    for j = i+1:n
        if size(E(E.i == i & E.j == j,:),1) > 1
            error('This table has multiedges')
        end
    end
end

Z = zeros(3*n);
A = zeros(n);
C = [];

% Identity constraints
for i = 1:n
    Z(3*i-2:3*i,3*i-2:3*i) = eye(3);
end

% Fill Z and A matrix
for k = 1:size(E,1)
    i = E(k,:).i;
    j = E(k,:).j;
    Z_ij = E(k,:).Z_ij;
    Z_ji = E(k,:).Z_ji;
    
    A(i,j) = 1;
    A(j,i) = 1;
    Z(3*i-2:3*i,3*j-2:3*j) = cell2mat(Z_ij);
    Z(3*j-2:3*j,3*i-2:3*i) = cell2mat(Z_ji);
end

% Fill C matrix
for k = 1:size(R,2)
    val = R{k};
    if size(val,1) > 0
        tic, pairs = nchoosek(val,2);
        for h = 1:size(pairs,1)
            a = pairs(h,1);
            b = pairs(h,2);
            
            c = zeros(1,3*n);
            c(1,3*a-2) = 1;
            c(1,3*b-2) = -1;
            C = [C;c];
            
            c = zeros(1,3*n);
            c(1,3*a-1) = 1;
            c(1,3*b-1) = -1;
            C = [C;c];
            
            c = zeros(1,3*n);
            c(1,3*a-0) = 1;
            c(1,3*b-0) = -1;
            C = [C;c];
        end
    end
end

C = C';
C = C / norm(C);
end