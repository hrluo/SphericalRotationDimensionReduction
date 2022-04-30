function distance_vec = MATCH_DIST(data_matrix,reduced_data_matrix,norm_type)
  %Function: MATCH_DIST
  %Hengrui Luo, Nov 07, 2020.
  %Functionality: Find the matched one-to-one distances (i.e., d(y_i,\hat{y}_i)) between dataset y and dimensionally reduced dataset \hat{y} in R^d (represented by n*d matrix data_matrix).
  %CAUTION: the distance is computed when both datasets are embedded in R^d, the original space with full dimension.
  %Input: data_matrix -- An n by d matrix representing n datapoints (rows) in R^d (columns)
  %       reduced_data_matrix -- An n by d matrix representing n datapoints (rows) in R^d (columns); d'<d but it should be a lower dimensional subspace embedded in R^d.
  %       norm_type -- options that indicates what kind of norms we use for calculating the pariwise distance, 
  %       it should be one of : 'L1', 'L2', 'Linf'
  %Output: distance_vec -- An 1 by n upper triangular matrix representing n matched pair distances between data_matrix and reduced_datamatrix
    s_raw = size(data_matrix);
    n_raw = s_raw(1);
    d_raw = s_raw(2);
    s_reduced = size(reduced_data_matrix);
    n_reduced = s_reduced(1);
    d_reduced = s_reduced(2);
    distance_vec = zeros(1,n_raw);
    if n_raw~=n_reduced
        fprintf('ERROR: The original dataset and the dimension reduced dataset must be of the same size.')
    end
    if d_raw~=d_reduced
        fprintf('ERROR: The original dataset and the reduced dataset should live in the same space R^d.')
        return
    end 
    for k = 1:n_raw
        if strcmpi(norm_type,'L1')
            distance_vec(k) = sum( abs( reduced_data_matrix(k,:) - data_matrix(k,:) ) );
        elseif strcmpi(norm_type,'L2')
            distance_vec(k) = sqrt( sum( ( reduced_data_matrix(k,:) - data_matrix(k,:) ).^2) );
        elseif strcmpi(norm_type,'Linf')
            distance_vec(k) = max( abs( reduced_data_matrix(k,:) - data_matrix(k,:) ) );
        else
        end
    end
    return
end