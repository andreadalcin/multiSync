function  x = inv_star(X)
if (nnz(X+X')~=0)
    X=(X-X')/2;
end

x=[X(3,2) X(1,3) X(2,1)]';

end
