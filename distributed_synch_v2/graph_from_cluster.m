function Z_c = graph_from_cluster(Z, A, map)
n = size(map,1);
Z_c = eye(3 * n);
for i = 1:n-1
    for j = i+1:n
        if A(i,j) ~= 0
            r_i = map(i);
            r_j = map(j);
            Z_c(3*i-2:3*i,3*j-2:3*j) = Z(3*r_i-2:3*r_i,3*r_j-2:3*r_j);
            Z_c(3*j-2:3*j,3*i-2:3*i) = Z(3*r_j-2:3*r_j,3*r_i-2:3*r_i);
        end
    end
end
end