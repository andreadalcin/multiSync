
function[Q,W,B,score] = Weighted_LAA(I,Q,QQ,Amatrix,Weights)
            N=max(max(I));
            w=zeros(size(QQ,1),4);W=zeros(N,4);
            i=I(1,:);j=I(2,:);

            % w=Qij*Qi
            w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
                repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
                [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

            % w=inv(Qj)*w=inv(Qj)*Qij*Qi
            w(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
                repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
                [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms


            s2=sqrt(sum(w(:,2:4).*w(:,2:4),2));
            w(:,1)=2*atan2(s2,w(:,1));
            i=w(:,1)<-pi;  w(i,1)=w(i,1)+2*pi;  i=w(:,1)>=pi;  w(i,1)=w(i,1)-2*pi;
            B=w(:,2:4).*repmat(w(:,1)./s2,[1,3]);
        

            B(isnan(B))=0;% This tackles the devide by zero problem.

            W(1,:)=[1 0 0 0];
           
            W(2:end,2:4)=(sparse(1:length(Weights),1:length(Weights),Weights,length(Weights),length(Weights))*Amatrix)\(repmat(Weights,[1,size(B,2)]).*B);
            
            score=sum(sqrt(sum(W(2:end,2:4).*W(2:end,2:4),2)))/N;

            theta=sqrt(sum(W(:,2:4).*W(:,2:4),2));
            W(:,1)=cos(theta/2);
            W(:,2:4)=W(:,2:4).*repmat(sin(theta/2)./theta,[1,3]);

            W(isnan(W))=0;

            Q=[ (Q(:,1).*W(:,1)-sum(Q(:,2:4).*W(:,2:4),2)),...  %scalar terms
                repmat(Q(:,1),[1,3]).*W(:,2:4) + repmat(W(:,1),[1,3]).*Q(:,2:4) + ...   %vector terms
                [Q(:,3).*W(:,4)-Q(:,4).*W(:,3),Q(:,4).*W(:,2)-Q(:,2).*W(:,4),Q(:,2).*W(:,3)-Q(:,3).*W(:,2)] ];   %cross product terms
            
end
