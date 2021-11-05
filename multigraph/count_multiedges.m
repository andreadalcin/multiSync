function count = count_multiedges(E, i, N)
count = 0;
E_i = E(E.i == i,:);

for j = i+1:N
    if size(E_i(E_i.j == j,:),1) > 1
        count = count + 1;
    end
end

end