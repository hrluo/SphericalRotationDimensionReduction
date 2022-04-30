clear
% For reproducibility
rng(123,'twister')  


T=readtable('BankNote.csv');
banknote = T{:,2:6};
%load('banknote.mat')
data = banknote;
X = banknote(:,1:4);
label = data(:,5);
label_set = unique(label);
% Set the sample size n and dimension d of the dataset for later use.
[n,d] = size(X);
% Set the target dimension.
retain_d = 2;

center1 = NaN;
radius1 = NaN;
W_matrix = eye(d);



%% Projected data

% SRCA
[output_SRCA,rotate_SRCA,opt_ind,center_SRCA,radius_SRCA,reduced_matrix] = SRCA(X,retain_d+1,'ALG',W_matrix,false,0,'PCA');


theta = zeros(n,1);
phi = zeros(n,1);
for i = 1:n
   phi(i) = real(atan2(reduced_matrix(i,2),reduced_matrix(i,1)));
   theta(i) = real(acos(reduced_matrix(i,3)/radius_SRCA));
%     theta(i) = theta(i)+0.8;
%     if theta(i) > pi/2
%        theta(i) = theta(i)-pi;
%     end
%     phi(i) = phi(i);
%     if phi(i) > pi
%        phi(i) = phi(i)-2*pi;
%     end
end
output_SRCA_intrinsic = [theta,phi];
plotClass(output_SRCA_intrinsic.',label)

% SPCA
[c,V,r]=Spherelets(X,retain_d);
output_SPCA = zeros(n,d);
output_SPCA_intrinsic = zeros(n,retain_d+1);
for i = 1:n
    output_SPCA_intrinsic(i,:) = (V.'*c+r*V.'*(X(i,:).'-c)/norm(V.'*(X(i,:).'-c))).';
    output_SPCA(i,:) = c.'+r*(X(i,:)-c.')*V*V.'/norm(V.'*(X(i,:).'-c));
end

%output_SPCA = output_SPCA_intrinsic;



theta = zeros(n,1);
phi = zeros(n,1);
for i = 1:n
    phi(i) = real(atan2(output_SPCA_intrinsic(i,2),output_SPCA_intrinsic(i,1)));
    theta(i) = real(acos(output_SPCA_intrinsic(i,3)/r));
    theta(i) = theta(i)+0.8;
    if theta(i) > pi/2
        theta(i) = theta(i)-pi;
    end
    phi(i) = phi(i)+5;
    if phi(i) > pi
       phi(i) = phi(i)-2*pi;
    end
end
output_SPCA_intrinsic = [theta,phi];
%plotClass(output_SPCA_intrinsic.',label)


% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(X);
output_PCA = score(:,1:retain_d);

% tSNE
output_tSNE = tsne(X);
%plot(output_tSNE(:,1),output_tSNE(:,2),'o')

% LLE
k = ceil(sqrt(n));
[output_LLE, mapping] = lle(X, retain_d, k);

% UMAP
[output_UMAP,umap,clusterIdentifiers,extras]=run_umap(X);



%% cluster preserving
figure


subplot(2,3,1)
plotClass(output_SRCA_intrinsic.',label)
title('SRCA')

subplot(2,3,2)
plotClass(output_SPCA_intrinsic.',label)
title('SPCA')

subplot(2,3,3)
plotClass(output_PCA.',label)
title('PCA')


subplot(2,3,4)
plotClass(output_tSNE.',label)
title('tSNE')

subplot(2,3,5)
plotClass(output_LLE.',label)
title('LLE')


subplot(2,3,6)
plotClass(output_UMAP.',label)
title('UMAP')

sgtitle('Cluster perserving for BankNote')


csvwrite('Banknote_output_PCA.csv', output_PCA)
csvwrite('Banknote_output_SRCA.csv', output_SRCA)
csvwrite('Banknote_output_SPCA.csv', output_SPCA)
csvwrite('Banknote_output_LLE.csv', output_LLE)
csvwrite('Banknote_output_tSNE.csv', output_tSNE)
csvwrite('Banknote_output_UMAP.csv', output_UMAP)

