function M = single_averaging_SO3_array(M_array)

% only one element: no averaging
if length(M_array) == 1
    M = M_array{1};
    return
end

% Array to matrix
M_mat = zeros(3,3,length(M_array));
for k = 1:length(M_array)
    M_mat(:,:,k) = M_array{k};
end

% Sum of the entries
M = mean(M_mat,3);

% Project onto SO(3)
[U,~,V]=svd(M);
M = U*diag([1 1 det(U*V')])*V';

end

