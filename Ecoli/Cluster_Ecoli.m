load('Ecoli.mat')
data = Ecoli;
X = data(:,3:9);
X = X{:,:};

label = data(:,10);
label = label{:,:};
label_set = unique(label);


label_small = [];
for j = 1:8
    if sum(label==label_set(j))<10
        label_small = [label_small,label_set(j)];
    end
end

for j = 1:size(label_small,2)
    X(find(label==label_small(j)),:)=[];
    label(find(label==label_small(j)))=[];
end


[n,d] = size(X);
retain_d = 2;

label_set_new = unique(label);
label_num = zeros(n,1);
for j=1:size(label_set_new,1)
    label_num(find(label==label_set_new(j)),1)=j*ones(sum(label==label_set_new(j)),1);
end
label = label_num;


center1 = NaN;
radius1 = NaN;
W_matrix = eye(d);


%% Projected data

% SPCA
[output_SPCA,rotate_SPCA,opt_ind,center_SPCA,radius_SPCA,reduced_matrix] = SRCA(X,retain_d+1,'ALG',W_matrix,false,0,'PCA');


theta = zeros(n,1);
phi = zeros(n,1);
for i = 1:n
   phi(i) = real(atan2(reduced_matrix(i,2),reduced_matrix(i,1)));
   theta(i) = real(acos(reduced_matrix(i,3)/radius_SPCA));
    phi(i) = phi(i)+1.7;
    if phi(i) > pi
       phi(i) = phi(i)-2*pi;
    end
end
output_SPCA_intrinsic = [theta,phi];
%plotClass(output_SPCA_intrinsic.',label)

% Spherelets
[c,V,r]=Spherelets(X,retain_d);
output_Spherelets = zeros(n,d);
output_Spherelets_intrinsic = zeros(n,retain_d+1);
for i = 1:n
    output_Spherelets_intrinsic(i,:) = (V.'*c+r*V.'*(X(i,:).'-c)/norm(V.'*(X(i,:).'-c))).';
    output_Spherelets(i,:) = c.'+r*(X(i,:)-c.')*V*V.'/norm(V.'*(X(i,:).'-c));
end

theta = zeros(n,1);
phi = zeros(n,1);
for i = 1:n
    phi(i) = real(atan2(output_Spherelets_intrinsic(i,2),output_Spherelets_intrinsic(i,1)));
    theta(i) = real(acos(output_Spherelets_intrinsic(i,3)/r));
    theta(i) = theta(i)+0.8;
    if theta(i) > pi/2
        theta(i) = theta(i)-pi;
    end
    phi(i) = phi(i)+5;
    if phi(i) > pi
       phi(i) = phi(i)-2*pi;
    end
end
output_Spherelets_intrinsic = [theta,phi];
plotClass(output_Spherelets_intrinsic.',label)


% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(X);
output_PCA = score(:,1:retain_d);



k = ceil(sqrt(n));


% tSNE
output_tSNE = tsne(X);



% LLE
[output_LLE, mapping] = lle(X, retain_d, k);


% UMAP
[output_UMAP,umap,clusterIdentifiers,extras]=run_umap(X);




%% cluster preserving
figure


subplot(2,3,1)
plotClass(output_SPCA_intrinsic.',label)
title('SPCA')

subplot(2,3,2)
plotClass(output_Spherelets_intrinsic.',label)
title('Spherelets')

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

sgtitle('Cluster perserving fo Ecoli')




