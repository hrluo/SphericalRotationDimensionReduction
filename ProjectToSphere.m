function [output_matrix,center,radius] = ProjectToSphere(data_matrix,center,radius,norm_type)
  %Function: Project the matrix to the same dimensional sphere described by
  %the given center and radius, you cannot project to a lower dimension.
  data_dim = size(data_matrix);
  output_matrix = 0*data_matrix;
  unit_vec = data_matrix;
  %use this matrix to store unit vector of a difference vector. 
  for k = 1:data_dim(1)
    unit_vec(k,:) = data_matrix(k,:) - center;
    unit_vec(k,:) = unit_vec(k,:)/norm(unit_vec(k,:),norm_type);
    output_matrix(k,:) = center + radius*unit_vec(k,:);
  end
end