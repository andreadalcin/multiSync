function EdgeTable = graph2edgetable(Z,A)
EdgeTable = empty_edgetable();
n = size(A,1);
for i = 1:n-1
    for j = i:n
        if A(i,j) ~= 0
            Z_ij = Z(3*i-2:3*i,3*j-2:3*j);
            Z_ji = Z(3*j-2:3*j,3*i-2:3*i);
            T = table(i,j,mat2cell(Z_ij,3,3),mat2cell(Z_ji,3,3),'VariableNames',{'i','j','Z_ij','Z_ji'});
            EdgeTable = [EdgeTable;T];
        end
    end
end
end