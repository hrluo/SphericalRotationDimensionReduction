clear
clc
% For reproducibility
rng(123,'twister')  
retain_d = 4;


T = importdata('climate.dat');
data = T(:,3:20);
Y = data;
[n,d] = size(Y);

%%%%%SRCA
W_matrix = eye(d);
global funcCallCount;
funcCallCount = 0; % Initialize the counter
[output_SRCA1,rotate_SRCA1,opt_ind1,center_SRCA1,radius_SRCA1,reduced_SRCA1] = hl_SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'PCA',1);
funcCallCount1 = funcCallCount;

funcCallCount = 0;
[output_SRCA2,rotate_SRCA2,opt_ind2,center_SRCA2,radius_SRCA2,reduced_SRCA2] = hl_SRCA(Y,retain_d+1,'ALG',W_matrix,true,0,'PCA',1);
funcCallCount2 = funcCallCount;

funcCallCount = 0;
[output_SRCA3,rotate_SRCA3,opt_ind3,center_SRCA3,radius_SRCA3,reduced_SRCA3] = SRCA(Y,retain_d+1,'ALG',W_matrix,false,0,'PCA');
funcCallCount3 = funcCallCount;

% funcCallCount = 0;
% [output_SRCA4,rotate_SRCA4,opt_ind4,center_SRCA4,radius_SRCA4,reduced_SRCA4] = SRCA(Y,retain_d+1,'ALG',W_matrix,true,0,'PCA');
% funcCallCount4 = funcCallCount;


disp(['(L2 binary search) number of evaluations = ',num2str(funcCallCount3)])
MD_SRCA3 = MATCH_DIST(Y,output_SRCA3,'L2').^2;
%disp(['opt_ind = ', num2str(opt_ind3),'; Estimated center:',num2str(center_SRCA3), '\t Estimated radius:',num2str(radius_SRCA3)])
display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA3))])


disp(['(L1 penalty)       number of evaluations = ',num2str(funcCallCount2)])
MD_SRCA2 = MATCH_DIST(Y,output_SRCA2,'L2').^2;
%disp(['opt_ind = ', num2str(opt_ind2),'; Estimated center:',num2str(center_SRCA2), '\t Estimated radius:',num2str(radius_SRCA2)])
display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA2))])


disp(['(branch-and-bound) number of evaluations = ',num2str(funcCallCount1)])
MD_SRCA1 = MATCH_DIST(Y,output_SRCA1,'L2').^2;
%disp(['opt_ind = ', num2str(opt_ind1),'; Estimated center:',num2str(center_SRCA1), '\t Estimated radius:',num2str(radius_SRCA1)])
display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA1))])



% disp(['(L1 penalty old)   number of evaluations = ',num2str(funcCallCount4)])
% MD_SRCA4 = MATCH_DIST(Y,output_SRCA4,'L2').^2;
% disp(['opt_ind = ', num2str(opt_ind4),'; Estimated center:',num2str(center_SRCA4), '\t Estimated radius:',num2str(radius_SRCA4)])
% display(['retain_dim = ',num2str(retain_d),' MSE of SRCA = ',num2str(mean(MD_SRCA4))])
% 
