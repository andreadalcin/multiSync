function [ Rm ] = dqq_rotation_quaternion_initialization( R )
QR=zeros(size(R,3),4);
for i=1:size(R,3)
    QR(i,:) = iquat(R(:,:,i));
end
SR=R(:,:,1);
SQ=iquat(SR);
for i=1:size(R,3)
    if norm(QR(i,:)+SQ)<norm(QR(i,:)-SQ)
        QR(i,:)=-QR(i,:);
    end
end
barQR=sum(QR,1)/size(R,3);
barQR=barQR/norm(barQR);
Rm=quat(barQR);
end






