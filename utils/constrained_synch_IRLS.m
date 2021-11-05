function R = constrained_synch_IRLS(Z,A,C,corr,maxiter,tol,cost_function,h)

n = size(A,1);
k = 1;
deltaW = 2 * tol;

W = A;
W_old = W;

while k <= maxiter && deltaW > tol
    %k
    R_array = constrained_synch_no_merge(Z,A,C,corr);
    R_mat = zeros(3,3,length(R_array));
    for i = 1:length(R_array)
        R_mat(:,:,i) = real(R_array{i});
    end
    
    W = update_weights_SO3(Z,R_mat,A,cost_function,h);
    
    deltaW = norm(W_old - W, 'fro') / (norm(W,'fro') * n);
    W_old = W;
    
    k=k+1;
end

% if use_bisquare
%     W = update_weights_SO3(Z,R,A,'bisquare',h,use_mex);
% end

% Merge replicas and correct for index differences
N = size(corr,2);

R = cell(1,N);
for i = 1:N
    if size(corr{i},2) > 0
        R{i} = R_array{corr{i}(:,1)};
    else
        R{i} = R_array{i};
    end
end

end

function [R] = constrained_synch_no_merge(Z,A,C,R)
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
R = cell(1,n);
for i = 1:n
    [U,~,V] = svd(Q(3*i-2:3*i,:));
    R{i} = U * diag([1,1,det(U*V')]) * V';
end
end

function [D] = update_weights_SO3(Z,R,A,weight_fun,h)

ncams=size(A,1);
[I,J]=find(triu(A,1));

R = reshape(permute(R,[1,3,2]),[],3);

T_omega = ((Z - (R * R') .* repelem(A,3,3)) .^2);
B0 = repelem(speye(ncams),1,3);
blknorm = sqrt(B0 * T_omega * B0');

[~,~,res]=find(triu(blknorm,1));

s =  mad(res(:),1)/0.6745;

weights = weightfun(res,s,weight_fun,h);

D = sparse([I;J;(1:ncams)'],[J;I;(1:ncams)'],[weights;weights;max(weights)*ones(ncams,1)],ncams,ncams);

end