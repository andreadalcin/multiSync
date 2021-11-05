function [R]=SO3_EIG(G,A)

ncams=size(A,1);

G = G.*repelem(A,3,3);

D = kron(diag(1./sum(A,2)),eye(3));

[M,~]=eigs(D*G,3); 

% Projection onto SO(3)
R=zeros(3,3,ncams);
for i=1:ncams
    [U,~,V] = svd(M(3*i-2:3*i,:)); 
    R(:,:,i)=U*V';
    if (det(R(:,:,i))<0)
        R(:,:,i)=-U*V';
    end
end


end