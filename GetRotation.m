function [rotate_matrix] = GetRotation(data_matrix,rotationMethod,retain_dim)
    %addpath('../src')
    data_dim = size(data_matrix);
    d = data_dim(2);
    fprintf('\n Rotation Method:', rotationMethod);
    
    if strcmpi(rotationMethod,'PCA')
        [coeff,score,latent] = pca(data_matrix,'Economy',false);
        rotate_matrix=coeff;%rotatefactors(coeff);
    end
    if strcmpi(rotationMethod,'ICA')
        Mdl = rica(data_matrix,d,'IterationLimit',100)
        rotate_matrix = Mdl.TransformWeights;
        if size(rotate_matrix,1)<size(rotate_matrix,2)
            rotate_matrix = transpose(rotate_matrix)
        end
    end
    if strcmpi(rotationMethod,'PCA_transpose')
        [coeff,score,latent] = pca(transpose(data_matrix));
        rotate_matrix=coeff;%rotatefactors(coeff);
    end
    if strcmpi(rotationMethod,'PCA_new')
        X_d = size(data_matrix);
        [V,lambda,mu] = pca_new(data_matrix,retain_dim+1);
        rotate_matrix=V*V.';
    end
    if strcmpi(rotationMethod,'RandPCA')
        data_dim=size(data_matrix);
        [coeff,score,latent] = rPCA(data_matrix,data_dim(2));
        rotate_matrix=coeff;
    end
    if strcmpi(rotationMethod,'FastPCA')
        data_dim=size(data_matrix);
        [coeff,score,latent] =  pcafast(data_matrix,data_dim(2));
        rotate_matrix=coeff;
    end
    if strcmpi(rotationMethod,'AdaptivePCA')
        data_dim=size(data_matrix);
        [coeff,score,latent] = adaptivepca(data_matrix,1.0d-13,data_dim(2));
        rotate_matrix=coeff;
    end
    if strcmpi(rotationMethod,'varimax')
        [X_factor_rotated, rotate_matrix] = rotatefactors(data_matrix,'Method','varimax','maxit',2000*d);
    end
    if strcmpi(rotationMethod,'orthomax')
        [X_factor_rotated, rotate_matrix] = rotatefactors(data_matrix,'Method','orthomax','maxit',2000*d);
    end
    if strcmpi(rotationMethod,'quartimax')
        [X_factor_rotated, rotate_matrix] = rotatefactors(data_matrix,'Method','quartimax','maxit',2000*d);
    end
    if strcmpi(rotationMethod,'equamax')
        [X_factor_rotated, rotate_matrix] = rotatefactors(data_matrix,'Method','equamax','maxit',2000*d);
    end
    if strcmpi(rotationMethod,'parsimax')
        [X_factor_rotated, rotate_matrix] = rotatefactors(data_matrix,'Method','parsimax','maxit',2000*d);
    end
    
    fprintf('\n Rotation Matrix Size:');
    size(rotate_matrix)
    
end
