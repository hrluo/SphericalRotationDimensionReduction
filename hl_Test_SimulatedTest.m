clear 
% For reproducibility
rng(123,'twister')  

% Set the sample size n and dimension d of the dataset for later use.
n = 1000;
d = 3;
% Set the target dimension.
retain_d = 2;

%Simulate a dataset on T2 of size n.
mu = zeros(1,d-1);
mu(1) = 1;
kappa = 0;
c_true = zeros(1,d);
r_true = 2;
R1=1/2;
R2=1/3;
Y = zeros(n,d);
for k = 1:n
    th = 2*pi*rand(1);
    ph = 2*pi*rand(1);
    Y(k,:) = [ (R1+R2*cos(th))*cos(ph), (R1+R2*cos(th))*sin(ph),R2*sin(th) ];     
    Y(k,:) = Y(k,:)* r_true;
    Y(k,:) = Y(k,:)+ c_true;
end

%%%%%SRCA
W_matrix = eye(d);
%[output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA, reduced_SRCA] = SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'PCA')
[output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA, reduced_SRCA, num_eval_SRCA] = hl_SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'PCA',1);
%[output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA, reduced_SRCA] = SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'PCA')
%%%%%Spherlets
[c_Spherelets,V,r_Spherelets]=Spherelets(Y,retain_d);
output_Spherelets = zeros(n,d);
for i = 1:n
    output_Spherelets(i,:) = c_Spherelets.'+ r_Spherelets*(Y(i,:)-c_Spherelets.')*V*V.'/norm((Y(i,:)-c_Spherelets.')*V*V.');
end
% PCA
[coeff,score,latent,tsquared,explained,mu_PCA] = pca(Y);
output_PCA = ones(n,1)*mu_PCA + (Y-ones(n,1)*mu_PCA)*coeff(:,1:retain_d)*coeff(:,1:retain_d).';

% plot3(0,0,0,'go')
% hold on
% plot3(Y(:,1),Y(:,2),Y(:,3),'r.')
% %plot3(Y_proj(:,1),Y_proj(:,2),Y_proj(:,3),'y.')
% plot3(output_SRCA(:,1),output_SRCA(:,2),output_SRCA(:,3),'k*')
% plot3(output_Spherelets(:,1),output_Spherelets(:,2),output_Spherelets(:,3),'bo')
% hold off
% 
% display(['true center is: ',num2str(c_true)])
% display(['Estimated center by Spherelets is: ',num2str(c_Spherelets.')])
% display(['Estimated center by SRCA is: ',num2str(center_SRCA)])
% display(['Number of binary evaluations used by SRCA is: ',num2str(num_eval_SRCA)])
% display(['true radius is ',num2str(r_true)])
% display(['Estimated radius by Spherelets is ',num2str(r_Spherelets)])
% display(['Estimated radius by SRCA is ',num2str(radius_SRCA)])
% 
% %matched pairwise distance distribution
% figure;
% hold on
% h_1 = histogram(MATCH_DIST(Y,output_Spherelets,'L2'),'BinWidth',0.01,'EdgeColor','none')
% h_2 = histogram(MATCH_DIST(Y,output_SRCA,'L2'),'BinWidth',0.01,'EdgeColor','none')
% h_3 = histogram(MATCH_DIST(Y,output_PCA,'L2'),'BinWidth',0.01,'EdgeColor','none')
% legend('Spherlets','SRCA','PCA')
% title('Matched pairwise distances of reduced-to-original histogram')
% hold off
% 
% MD_PCA = MATCH_DIST(Y,output_PCA,'L2').^2;
% MD_SRCA = MATCH_DIST(Y,output_SRCA,'L2').^2;
% MD_Spherelets= MATCH_DIST(Y,output_Spherelets,'L2').^2;
% display(['retain_dim = ',num2str(retain_d),' MSE of PCA = ',num2str(mean(MD_PCA))])
% display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA))])
% display(['retain_dim = ',num2str(retain_d),' MSE of Spherelets = ',num2str(mean(MD_Spherelets))])