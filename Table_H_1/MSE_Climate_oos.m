clear
% For reproducibility
rng(123,'twister')

T = importdata('climate.dat');
data = T(:,3:20);
Y = data;
% Set the sample size n and dimension d of the dataset for later use.
[n,d] = size(Y);

W_matrix = eye(d);

MSE_PCA = zeros(1,d-1);
MSE_SRCA = zeros(1,d-1);
MSE_SPCA = zeros(1,d-1);

idx = randperm(n);
ntr = ceil(n/2);
nte = n - ntr;
Ytr = Y(idx(1:ntr),:);
Yte = Y(idx(ntr+1:n),:);

for retain_d = 1:4
    
    %SRCA
    %tic
    [output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA,reduced_SRCA] = SRCA(Ytr,retain_d+1,'ALG',W_matrix,true,1,'ica');
    %toc
    Proj_SRCA = SRCA_projection(Yte,center_SRCA,radius_SRCA,rotate_SRCA,opt_ind);
    
    
    
    
    %%%%%SPCA
    %tic
    [c_SPCA,V,r_SPCA]=Spherelets(Ytr,retain_d);
    %toc    
    Proj_SPCA = zeros(nte,d);
    for i = 1:nte
        Proj_SPCA(i,:) = c_SPCA.'+ r_SPCA*(Yte(i,:)-c_SPCA.')*V*V.'/norm((Yte(i,:)-c_SPCA.')*V*V.');
    end
    
    % PCA
    [coeff,score,latent,tsquared,explained,mu_PCA] = pca(Ytr);
    Proj_PCA = ones(nte,1)*mu_PCA + (Yte-ones(nte,1)*mu_PCA)*coeff(:,1:retain_d)*coeff(:,1:retain_d).';
    
    
    
    MD_PCA = MATCH_DIST(Yte,Proj_PCA,'L2').^2;
    MD_SPCA= MATCH_DIST(Yte,Proj_SPCA,'L2').^2;
    MD_SRCA = MATCH_DIST(Yte,Proj_SRCA,'L2').^2;
    
    MSE_PCA(retain_d) = mean(MD_PCA);
    MSE_SPCA(retain_d) = mean(MD_SPCA);
    MSE_SRCA(retain_d) = mean(MD_SRCA);
    
    display(['retain_dim = ',num2str(retain_d),' MSE of PCA = ',num2str(mean(MD_PCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SPCA = ',num2str(mean(MD_SPCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA))])
    
end

MSEs = [MSE_PCA; MSE_SPCA; MSE_SRCA]

% figure
% hold on
% plot(1:1:d-1,MSE_PCA,'r','LineWidth',4)
% plot(1:1:d-1,MSE_SRCA,'k','LineWidth',4)
% plot(1:1:d-1,MSE_SPCA,'b','LineWidth',4)
% hold off
% legend('PCA','SRCA','SPCA','FontSize',100)
% xlabel('d''','FontSize',100)
% ylabel('MSE','FontSize',100)
% title('MSEs of Climate','FontSize',100)





