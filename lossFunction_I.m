function LOSS = lossFunction_I(data_matrix,ind_vec,cov_mat,lambda,center,radius)
    LOSS = 0;
    data_dim = size(data_matrix);
    id_I = diag(ind_vec);
    %disp(cov_mat)
    %disp(size(id_I))
    %%disp(data_dim)
    %id_Ic = eye(data_dim(2))-id_I;
    %for i = 1:data_dim(1)
    %  X_i = data_matrix(i,:);
    %  %LOSS = LOSS + (X_i - center)*transpose(X_i - center) + radius*radius - 2*radius*sqrt( (X_i - center)*cov_mat*id_I*transpose(X_i - center) ) + lambda*sum(abs(ind_vec));
    %  LOSS = LOSS + (X_i - center)*transpose(X_i - center) + radius*radius  - 2*radius*sqrt( (X_i - center)*cov_mat*id_I*transpose(X_i - center) );
    %  %Didong's formulation.
    %  %LOSS = LOSS + (X_i - center)*id_I*transpose(X_i - center) + (sqrt( (X_i - center)*id_Ic*transpose(X_i - center) )-radius)^2 + lambda*sum(abs(ind_vec));
    %  %Original formulation.
    %end
    
    LOSS = sum( diag( (data_matrix - center)*cov_mat*transpose(data_matrix - center) ) );
    %disp('dimensionss>')
    %size(cov_mat)
    %size(data_matrix)
    %size(center)
    LOSS = LOSS + radius*radius*data_dim(1);
    LOSS = LOSS -  2*radius*sum( sqrt( diag( (data_matrix - center)*cov_mat*id_I*transpose(data_matrix - center) ) ) );
   
    LOSS = LOSS + lambda*sum( abs(ind_vec) );  
    return;