function E = sort_edges(E)
for k = 1:size(E,1)
    E_k = E(k,:);
    i = E_k.i;
    j = E_k.j;
    if i > j
        Z_ij = E_k.Z_ij;
        Z_ji = E_k.Z_ji;
        E(k,:) = [j,i,Z_ji,Z_ij];
    end
end
tic
E = sortrows(E,{'i','j'});
end