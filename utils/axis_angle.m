function [theta,u,v] = axis_angle(R)
theta = acos((trace(R)-1)/2);
theta=real(theta); 
u = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]'; 
if (nnz(R-R')==0)
    if (nnz(R-eye(3))==0)
        theta=0;
        u=[1;1;1];
    else
        theta=pi;
        A=R+eye(3);
        u=A(:,1);
    end       
end
u=u/norm(u);
v=u*theta;
end
