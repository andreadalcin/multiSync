function  S = star(x)

if size(x,2) ~=1 x = x'; end

if (size(x,2) ~=1)
    error('Argument must be a vector');
end

n = length(x);
S =  [];

for i=1:n
    S = [S;
        zeros(n-i, i-1), -x(i+1:end), x(i) * eye(n-i)];
end

if n == 3
    t = S(3,:); S(3,:)=S(1,:); S(1,:)=t;
    S(2,:)=-S(2,:);
end


