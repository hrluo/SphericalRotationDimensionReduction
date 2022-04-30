clear
% For reproducibility
rng(123,'twister')
T=readtable('BankNote.csv');
banknote = T{:,2:6};
data = banknote;
Y = banknote(:,1:4);

[n,d] = size(Y);
W_matrix = eye(d);

MSE_PCA = zeros(1,d-1);
MSE_SPCA = zeros(1,d-1);
MSE_SRCA = zeros(1,d-1);
for retain_d = 1:d-1
    
    %SRCA
    %tic
    [output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA,reduced_SRCA] = SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'pca');
    %toc
    Proj_SRCA = SRCA_projection(Y,center_SRCA,radius_SRCA,rotate_SRCA,opt_ind);
    
    
    
    
    %%%%%SPCA
    %tic
    [c_SPCA,V,r_SPCA]=Spherelets(Y,retain_d);
    %toc    
    Proj_SPCA = zeros(n,d);
    for i = 1:n
        Proj_SPCA(i,:) = c_SPCA.'+ r_SPCA*(Y(i,:)-c_SPCA.')*V*V.'/norm((Y(i,:)-c_SPCA.')*V*V.');
    end
    
    % PCA
    [coeff,score,latent,tsquared,explained,mu_PCA] = pca(Y);
    Proj_PCA = ones(n,1)*mu_PCA + (Y-ones(n,1)*mu_PCA)*coeff(:,1:retain_d)*coeff(:,1:retain_d).';
    
    
    
    MD_PCA = MATCH_DIST(Y,Proj_PCA,'L2').^2;
    MD_SPCA= MATCH_DIST(Y,Proj_SPCA,'L2').^2;
    MD_SRCA = MATCH_DIST(Y,Proj_SRCA,'L2').^2;
    
    MSE_PCA(retain_d) = mean(MD_PCA);
    MSE_SPCA(retain_d) = mean(MD_SPCA);
    MSE_SRCA(retain_d) = mean(MD_SRCA);
    
    display(['retain_dim = ',num2str(retain_d),' MSE of PCA = ',num2str(mean(MD_PCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SPCA = ',num2str(mean(MD_SPCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA))])
    
end

MSEs = [MSE_PCA; MSE_SPCA; MSE_SRCA];
