function [R,k,W] = rotation_synch_IRLS(X,A,maxiter,tol,cost_function,h)

n = size(A,1);

k=1;
deltaW=2*tol;

W=A;
W_old = W;

while k<=maxiter && deltaW>tol
    %k
    [R_mat] = SO3_EIG(X,W);
    
    
    W=update_weights_SO3(X,R_mat,A,cost_function,h);
    
    deltaW = norm(W_old -W,'fro')/(norm(W,'fro')*n);
    W_old = W;
    
    k=k+1;
    
end

R = cell(1,n);
for i = 1:length(R)
    R{i} = R_mat(:,:,i);
end

end



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





function [D]=update_weights_SO3(X,R,A,weight_fun,h)

ncams=size(A,1);
[I,J]=find(triu(A,1));


R=reshape(permute(R,[1,3,2]),[],3);

T_omega = ( (X - (R*R').*repelem(A,3,3)).^2);
B0 = repelem(speye(ncams),1,3);
blknorm = sqrt(B0*T_omega*B0');

[~,~,res]=find(triu(blknorm,1));

s =  mad(res(:),1)/0.6745;

weights = weightfun(res,s,weight_fun,h);

D=sparse([I;J;(1:ncams)'],[J;I;(1:ncams)'],[weights;weights;max(weights)*ones(ncams,1)],ncams,ncams);

end