function [R,k,W] = SO3_eig_IRLS(X,A,maxiter,tol,cost_function,h,use_mex,use_bisquare)

if nargin<6
    use_mex=0;
end

n = size(A,1);

k=1;
deltaW=2*tol;

W=A;
W_old = W;

while k<=maxiter && deltaW>tol
    
    %k
    
    [R]=SO3_EIG(X,W);
    
    W=update_weights_SO3(X,R,A,cost_function,h,use_mex);
    
    deltaW = norm(W_old -W,'fro')/(norm(W,'fro')*n);
    W_old = W;
    
    k=k+1;
    
end

if use_bisquare
    W=update_weights_SO3(X,R,A,'bisquare',h,use_mex);
end

end



function [D]=update_weights_SO3(X,R,A,weight_fun,h,use_mex)

ncams=size(A,1);
[I,J]=find(triu(A,1));


if use_mex
    [res]=residuals_EIG_SO3_mex(I,J,full(X),R);
else
    R=reshape(permute(R,[1,3,2]),[],3);
    
    T_omega = ( (X - (R*R').*repelem(A,3,3)).^2);
    B0 = repelem(speye(ncams),1,3);   
    blknorm = sqrt(B0*T_omega*B0');
    
    [~,~,res]=find(triu(blknorm,1));

end

%s =  0.6745*mad(res(:),1);
s =  mad(res(:),1)/0.6745;

weights = weightfun(res,s,weight_fun,h);

D=sparse([I;J;(1:ncams)'],[J;I;(1:ncams)'],[weights;weights;max(weights)*ones(ncams,1)],ncams,ncams);

end