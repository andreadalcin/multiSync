function [G,A,R_gt,cameras_gt,i_cc,i_gt] = load_dataset(opts)

data_name = opts.data_name;
debug = opts.debug;

f=fopen(strcat('./datasets_snavely/',data_name,'/gt_bundle.out'));

a=fscanf(f,'%c',[18 1]);
disp(['Ground Truth: ',a'])
b=fscanf(f,'%u',[2 1]);
n=b(1);
ground_truth=fscanf(f,'%f',[3,5*n]);
ground_truth=ground_truth';

fclose(f);

R_gt=nan(3,3,n); 

cameras_gt=[]; 
for i=1:n
    
    Ri=ground_truth(5*i-3:5*i-1,:);
    Ti=ground_truth(5*i,:)';
    
    if all(all(Ri==0))
        R_gt(:,:,i)=NaN(3,3);
        T_gt(:,i)=NaN(3,1);
    else
        R_gt(:,:,i)=Ri;
        T_gt(:,i)=Ti;
        cameras_gt=[cameras_gt i-1];
    end
    
end


f=fopen(strcat('./datasets_snavely/',data_name,'/cc.txt'));
cc=fscanf(f,'%u',[1,inf]);
fclose(f);

cc=sort(cc);
ncams=length(cc);


f=fopen(strcat('./datasets_snavely/',data_name,'/EGs.txt'));

[EGs]=fscanf(f,'%f',[14 inf]);
EGs=EGs';
npairs=size(EGs,1);
fclose(f);

G=eye(3*ncams); 
A=eye(ncams); 
E_R=nan(ncams); 
R_gt_reduced=R_gt(:,:,cc+1); 

for p=1:npairs
    
    cam_i=EGs(p,1);
    cam_j=EGs(p,2);
    
    [ans_i,i]=find(cc==cam_i);
    [ans_j,j]=find(cc==cam_j);
    
    if ans_i && ans_j
        
        Rij=[EGs(p,3) EGs(p,4) EGs(p,5);EGs(p,6) EGs(p,7) EGs(p,8);EGs(p,9) EGs(p,10) EGs(p,11)];       
        
        G(3*i-2:3*i,3*j-2:3*j)=Rij;
        G(3*j-2:3*j,3*i-2:3*i)=Rij';
        
        A(i,j)=1;
        A(j,i)=1;
        
        if isfinite(R_gt_reduced(1,1,i)) &&  isfinite(R_gt_reduced(1,1,j))
            E_R(i,j)=phi_6(R_gt_reduced(:,:,i)*R_gt_reduced(:,:,j)',Rij)*180/pi;
        end
        
    end
end

G=sparse(G);
A=sparse(A);

fraction_missing=(1-(nnz(A)-ncams)/nchoosek(ncams,2)/2)*100;
disp(['Fraction of missing pairs = ',num2str(fraction_missing),' %'])

disp('Error [degrees] on pairwise rotations:')
disp(['Mean error = ',num2str(nanmean(E_R(:)))])
disp(['Median error = ',num2str(nanmedian(E_R(:)))])

disp(['Number of cameras = ',num2str(ncams)])

figure
hist(E_R(:),20)
title('Input (pairwise) rotations')
xlabel('Error [degrees]')
set(gca,'FontSize',16,'LineWidth',3)
grid on

[index,i_cc,i_gt]=intersect(cc,cameras_gt);

[~,ind]=max(sum(A));

D=diag(sum(A,2));
lambda2=max(abs(eigs(D-A,2,'sm')));
index_snavely=lambda2/ncams

end