function [X, Z, A] = random_graph_SLD(n,d)
X = zeros(d*n,d);
Z = eye(d*n);
A = zeros(n);

% Generate nodes in the graph
for i = 1:n
    Xi = rand(d);
    while det(Xi) <= 0
        Xi = rand(d);
    end
    Xi = Xi ./ nthroot(det(Xi),d);
    X(d*(i-1)+1:d*i,:) = Xi;
end
X = X / X(1:d,:);

% Compute consistency constraints
for i = 1:n-1
    for j = i+1:n
        A(i,j) = 1;
        A(j,i) = 1;
        
        X_i = X(d*(i-1)+1:d*i,:);
        X_j = X(d*(j-1)+1:d*j,:);
        
        Z(d*(i-1)+1:d*i,d*(j-1)+1:d*j) = X_i / X_j;
        Z(d*(j-1)+1:d*j,d*(i-1)+1:d*i) = X_j / X_i;
    end
end

end