function [ Rmean ] = L1_single_averaging( R,iter_max)

S(:,:,1) = dqq_rotation_quaternion_initialization( R );
nofR=size(R);

iter=1;


while isreal(S(:,:,iter)) && iter<=iter_max;
    
    iter=iter+1;

    sum_vmatrix_normed(:,:,iter)=zeros(3,3);
       
    conta=0;
    for j=1:size(R,3)
        
        [theta,u] = axis_angle(R(:,:,j)*S(:,:,iter-1)');

        vmatrix(:,:,j)=star(theta*u);
        
        vmatrix_normed(:,:,j)=vmatrix(:,:,j)/norm(vmatrix(:,:,j));
        sum_vmatrix_normed(:,:,iter)=sum_vmatrix_normed(:,:,iter)+vmatrix_normed(:,:,j);
        inv_norm_vmatrix(j)=1/norm(vmatrix(:,:,j));
        
        
        
        if norm(vmatrix(:,:,j))<1e-6
            conta=conta+1;
        end
        
    end
    


    
    if conta==size(R,3)
        S(:,:,iter)=S(:,:,iter-1);
        break
    end
    
    
    while conta>0 && conta<size(R,3)
        
        S(:,:,iter-1)=S(:,:,iter-1)+rand(3)*(1e-4);
        [u,~,v]=svd(S(:,:,iter-1));
        S(:,:,iter-1)=u*diag([1,1,det(u*v')])*v';
        
        conta=0;
        sum_vmatrix_normed(:,:,iter)=zeros(3,3);
        
        for j=1:size(R,3)
            
            [theta,u] = axis_angle(R(:,:,j)*S(:,:,iter-1)');
            vmatrix(:,:,j)=star(theta*u);
            
            
            vmatrix_normed(:,:,j)=vmatrix(:,:,j)/norm(vmatrix(:,:,j));
            sum_vmatrix_normed(:,:,iter)=sum_vmatrix_normed(:,:,iter)+vmatrix_normed(:,:,j);
            inv_norm_vmatrix(j)=1/norm(vmatrix(:,:,j));
            
            if norm(vmatrix(:,:,j))<1e-6
                conta=conta+1;
            end
        end
        
        
    end
    
    
    if conta==size(R,3)
        S(:,:,iter)=S(:,:,iter-1);
        break
    end
    
    delta(:,:,iter)=sum_vmatrix_normed(:,:,iter)/sum(inv_norm_vmatrix);
    
    
    u=inv_star(delta(:,:,iter));
    theta=norm(u);
    u=u/theta;
    S(:,:,iter) = inv_axis_angle(theta,u)*S(:,:,iter-1);
       
end

Rmean=S(:,:,iter);

end

