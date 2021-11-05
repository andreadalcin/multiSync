
function d=phi_6(R1,R2)

R=R1*R2';
d = acos((trace(R)-1)/2);

d=real(d);

end
