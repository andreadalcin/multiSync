function R = rotation_synch(Z,A)
n = size(A,1);
D = diag(1 ./ sum(A,2));

[Q,~] = eigs(kron(D, eye(3)) * Z, 3);

% Remove ambiguity
Q = Q / (Q(1:3,:));

% Projection onto SO(3)
R = cell(1,n);
% R = zeros(3*N,3);
for i = 1:n
    [U,~,V] = svd(Q(3*i-2:3*i,:));
    % R(3*i-2:3*i,:) = U * diag([1,1,det(U*V')]) * V';
    R{i} = U * diag([1,1,det(U*V')]) * V';
end
end