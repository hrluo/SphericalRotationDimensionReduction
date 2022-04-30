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

% Set the target dimension.
for retain_d = 1:d-1
    
    %%%%%SRCA
    %tic
    [output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA,reduced_SRCA] = SRCA(Y,retain_d+1,'ALG',W_matrix,true,0.1,'PCA');
    %toc
    Proj_SRCA = zeros(n,d);
  
    for i = 1:n
        Proj_SRCA(i,:) = center_SRCA+ radius_SRCA*(Y(i,:)-center_SRCA)*rotate_SRCA.*opt_ind*rotate_SRCA.'/norm((Y(i,:)-center_SRCA)*rotate_SRCA.*opt_ind*rotate_SRCA.');
    end
    
    
    
    %%%%%Spherlets
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
    MD_SRCA = MATCH_DIST(Y,Proj_SRCA,'L2').^2;
    MD_SPCA= MATCH_DIST(Y,Proj_SPCA,'L2').^2;
    
    MSE_PCA(retain_d) = mean(MD_PCA);
    MSE_SRCA(retain_d) = mean(MD_SRCA);
    MSE_SPCA(retain_d) = mean(MD_SPCA);
    
    display(['retain_dim = ',num2str(retain_d),' MSE of PCA = ',num2str(mean(MD_PCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA))])
    display(['retain_dim = ',num2str(retain_d),' MSE of SPCA = ',num2str(mean(MD_SPCA))])
    
end

MSEs = [MSE_PCA;MSE_SRCA;MSE_SPCA];

figure
hold on
plot(1:1:d-1,MSE_PCA,'r','LineWidth',4)
plot(1:1:d-1,MSE_SRCA,'k','LineWidth',4)
plot(1:1:d-1,MSE_SPCA,'b','LineWidth',4)
hold off
legend('PCA','SRCA','SPCA','FontSize',100)
xlabel('d''','FontSize',100)
ylabel('MSE','FontSize',100)
title('MSEs of Climate','FontSize',100)





