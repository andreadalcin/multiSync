function M = single_averaging_SO3(M_values)

% If only one element exists, no averaging
tic
if size(M_values,3) == 1
    M = M_values;
    return
end

M_values(isnan(M_values))=0;

% Sum of the entries
M = mean(M_values,3);

% Project onto SO(3)
[U,~,V] = svd(M);
M = U * diag([1 1 det(U*V')]) * V';
%toc

end