function [X, Z, A] = random_graph(N)
X = zeros(3*N,3);
Z = zeros(3*N);
A = zeros(N);

% Generate nodes in the graph
for n = 1:N
    rot = rot_y(rand * 2 * pi) * rot_x(rand * 2 * pi);
    X(3*n-2:3*n,:) = rot;
end
X = X / X(1:3,:);

% Identity constraints on the diagonal
for i = 1:N
    Z(3*i-2:3*i,3*i-2:3*i) = eye(3);
end

% Compute consistency constraints
for i = 1:N-1
    for j = i+1:N
        A(i,j) = 1;
        A(j,i) = 1;
        
        X_i = X(3*i-2:3*i,:);
        X_j = X(3*j-2:3*j,:);
        
        Z_ij = X_i / X_j;
        Z_ji = X_j / X_i;
        
        Z(3*i-2:3*i,3*j-2:3*j) = Z_ij;
        Z(3*j-2:3*j,3*i-2:3*i) = Z_ji;
    end
end

end