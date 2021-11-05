function R = inv_axis_angle(theta,u)
if (length(u) ~= 3)
    error('error');
end

if theta ==0
    u=[1;1;1];
end

u=u/norm(u); 
u=u(:);

R=cos(theta)*eye(3)+sin(theta)*star(u)+(1-cos(theta))*u*u';

end


