function synch_multi = constrained_synch(Z,A,C,R)
n = size(A,1);
D = kron(diag(sum(A,2)), eye(3));
B = Z - D;
P = eye(size(B,1)) - C * pinv(C);

[Q,O] = eigs(P*(B'*B),rank(C)+3,'smallestabs');

% Filter the smallest rank(C) eigenvalues
[~,ord] = sort(diag(O));
Q = Q(:,ord);
Q = Q(:,rank(C)+1:end);

% Remove ambiguity
idx = R{1}(1);
Q = Q / (Q(idx*3-2:idx*3,:));

% Projection onto SO(3)
synch = cell(1,n);
for i = 1:n
    [U,~,V] = svd(Q(3*i-2:3*i,:));
    synch{i} = U * diag([1,1,det(U*V')]) * V';
end

% Merge replicas and correct for index differences
N = size(R,2);

synch_multi = cell(1,N);
for i = 1:N
    if size(R{i},2) > 0
        synch_multi{i} = synch{R{i}(:,1)};
    else
        synch_multi{i} = synch{i};
    end
end

end