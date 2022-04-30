function [output_matrix,rotate_matrix,opt_ind,center,radius,reduced_matrix] = SRCA(data_matrix,retain_dim,loss_type,W_matrix,enable_norm,lambda,rotationMethod,optim_steps)
  %Function: SRCA-joint, based on Didong's suggestion of one-step method.
  %Hengrui Luo and Didong Li, Aug 1st, 2021.
  %This scipt implements to SRCA method.
  %Functionality: Find the optimal lower-dimensional sphere (specified by retain_dim) that represents the dataset in R^d (represented by n*d matrix data_matrix).
  %Input: data_matrix -- An n by d matrix representing n datapoints (rows) in R^d (columns)
  %       retain_dim -- It must be less or equal to d, the dimension of the reduced dataset.
  %       center -- An 1 by d vector. By default it is NULL, then we would estimate the center. If provided, we will use it to standardize our dataset.
  %       radius -- A scalar. By default it is NULL, then we would estimate the radius. If provided, we will use it to standardize our dataset.
  %       loss_type -- either 'ALG' for algebraic loss function for estimation or 'GEO' for geometric loss function. However, this option is useful only if center or radius is not provided.
  %       VERBOSE -- if TRUE, more information for output will be printed
  %       W_matrix -- By default it is NULL, meaning the identity matrix of dimension d would be used.
  %       enable_norm -- if true, we use L1 norm to approximate the binary search problem.
  %       lambda -- the penalty parameter lambda>0 to adjust result for the sparsity in the dataset.
  %       rotationMethod -- which rotation method you want to use for the
  %       preprocessing the data matrix. PCA/varimax/orthomax/quartimax/equamax/parsimax
  %
  % Full objective function for binary search: (expression for each single term x_i)
  %\begin{align*}
  %\min_{\bm{v}\in\mathbb{S}_{l_{k}}^{p}}\sum_{i=1}^{n}-2\sqrt{\bm{x}_{i}^{T}\sqrt{\bm{W}}^{T}\bm{v}^{T}\bm{I}_{p}\bm{v}\sqrt{\bm{W}}\bm{x}_{i}} & +\lambda\|\bm{I}_{\mathcal{I}}\bm{x}_{i}\|_{1},\end{align*}
  %\text{ s.t. }\|\bm{v}\|_{l_{k}}\leq p_{0},\lambda>0,
  %\end{align*}
  %Output: output_matrix -- An n by d (but it represents a retain_dim dimensional subspce) matrix representing n datapoints (rows) in R^d (columns)
  %        rotate_matrix -- An d by d matrix we use in SRCA, provided by PCA procedure.
  %        opt_ind -- A binary vector with 0 or 1 as its entries, representing the optimal solution to the binary optimization problem.
  %        center -- the center (estimate, if center=NULL in input) used in SRCA.
  %        radius -- the radius (estimate, if center=NULL in input) used in SRCA.
  %        reduced_matrix -- And n by retain_dim matrix matrix representing n datapoint centered at the origin. 
  %        
  % If you want to debug or have a closer look at the procedure, use
  % VERBOSE = true;
  %Testing varimax rotation instead of PCA rotation
    %rotationMethod = 'varimax';
    fprintf('SRCA-joint Method: Only algebraic loss function is implemented, ALG and GEO parameter would be ignored. \n')
    if nargin<7
        %By default we use binary search
        fprintf('SRCA-joint Method: Binary search...\n')
        L1_NORM = false;
    else
        L1_NORM = enable_norm;
        if L1_NORM == true 
            fprintf('SRCA-joint Method: L1 optimization...\n')
        else
            fprintf('SRCA-joint Method: Binary search...\n')
        end
    end    
    %Do we want more outputs for debug purpose?
    VERBOSE = true;
    %How many restarts shall we try during the optimization?
    RESTART = 1;
    %dimension of data matrix
    data_dim = size(data_matrix);
    fprintf('Data matrix dimensions:\n')
    disp(data_dim)
    d = data_dim(2);
    
    if retain_dim>data_dim(2)
        error('ERROR: The retained dimension cannot be greater than the dimension of the original data matrix.')
    end
    if nargin<8
        optim_steps = 500000;
    end
    if numel(sum(isnan(W_matrix))>0) && numel(L1_NORM == false)
        W_matrix=eye(data_dim(2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Step 1: Get the empirical mean for PCA .      %
    empirical_mean = mean(data_matrix,1);
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Step 2:Conduct the rotation given by PCA         % 
    rotate_matrix = GetRotation(data_matrix,rotationMethod,retain_dim);
    %disp(size(data_matrix-empirical_mean))
    %disp(size(rotate_matrix))
    X_rotated = (data_matrix-empirical_mean)*rotate_matrix;
    fprintf('Rotated dimensions:\n')
    d_rotated = size(X_rotated);
    disp(d_rotated)
    d_rotated = d_rotated(2);
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Step3: Binary search for best axes (dimension reduction)  %
    % Currently, this step only searches the optimal indices \mathcal{I}.
    if L1_NORM == false
        fprintf('SRCA-joint Method: Binary search start.\n')
        bin_list = dec2bin(2^d-1:-1:0)-'0';
        min_val = Inf;
        opt_ind = bin_list(1,:);
        for k = 1:2^d
            if sum(bin_list(k,:))~=retain_dim
                continue
            end
            disp(k)
            %Estimate center and radius jointly
            %lossfun = @(x)ALG_LOSS(x,data_matrix,W_matrix);
            %lossfun = @(x)GEO_LOSS(x,data_matrix,W_matrix);
            lossfun = @(x)lossFunction_I(X_rotated,bin_list(k,:),W_matrix,lambda,x(1:end-1),x(end));
            [center_candidate,radius_candidate] = EstimateSphere(X_rotated,lossfun);
            cur_val = lossFunction_I(X_rotated,bin_list(k,:),W_matrix,lambda,center_candidate,radius_candidate);
            if cur_val<min_val
                min_val = cur_val;
                cur_center = center_candidate;
                cur_radius = radius_candidate;
                opt_ind = bin_list(k,:);
            end
        end
        re_estimate=false;
    else
        fprintf('SRCA-joint Method: L1 optimization start.\n')
        fprintf('SRCA-joint Method: L1 only supports identity W matrix at the moment.\n')
        %if sum(isnan(W_matrix))>0
        %W_matrix=eye(d_rotated);
        %size(W_matrix)
        %end
        binary_obj_fun = @(x)lossFunction_I(X_rotated,x(1:d_rotated),W_matrix,lambda,x(d_rotated+1:end-1),x(end));
        lb = horzcat( ones(1,d_rotated)*0., ones(1,d_rotated)*(-inf),0 );
        ub = horzcat( ones(1,d_rotated)*1., ones(1,d_rotated)*(inf),inf);
        %lb = [];
        %ub = [];
        x0 = zeros(1,2*d_rotated+1);
        x0(end) = 1;
        A = ones(1,length(x0));%[];
        b = retain_dim;
        Aeq = [];
        beq = [];
        nonlcon =  @(x)nonlcon_init(x,retain_dim,d_rotated);
        %nonlcon =  @(x)sum(abs(x))-retain_dim;
        %The nonlinear constraint is nonlcon(x)\leq 0.
        options_con = optimoptions('fmincon','Algorithm','interior-point',...
                                    'HessianApproximation','bfgs',...
                                    'Display','iter-detailed',...%'final',...%
                                    'MaxFunctionEvaluations',optim_steps,...
                                    'FiniteDifferenceType','central',...
                                    'MaxIterations',2000,...
                                    'ScaleProblem',true,...
                                    'SubproblemAlgorithm','cg',...
                                    'ConstraintTolerance',1e-12,...
                                    'PlotFcn',[],...%'optimplotx' ,...
                                    'StepTolerance',1e-12,...
                                    'OptimalityTolerance',1e-12);
        [x_t,fval_t] = fmincon(binary_obj_fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options_con);

        opt_ind = x_t(1:d_rotated);
        %cur_center = x_t(d+1:end-1);
        %cur_radius = x_t(end);
        re_estimate = true;
    end
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
                                                              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %         Step 4: Eliminate the un-chosen dimensions       %
    [~,chosen_index]=maxk(abs(opt_ind), retain_dim);
    %Keep the k-largest coordinate coefficients.
    %Or we can set a threshold 0.01*1/d can be increase or reduced; especially in L1-enabled mode.
    chosen_mult = zeros(1,d_rotated);
    chosen_mult(chosen_index) = 1;
    %fprintf('chosen indices')
    %chosen_mult
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %         Step 5: Re-estimate the center and radius        %
    if re_estimate == true
        lossfun = @(x)lossFunction_I(X_rotated,chosen_mult,W_matrix,lambda,x(1:end-1),x(end));
        [center_candidate,radius_candidate] = EstimateSphere(X_rotated,lossfun);
        cur_center = center_candidate;
        cur_radius = radius_candidate;
    end
    if d==d_rotated
        %Check if the rotation leads to rank-deficient matrix, if so, do
        %not estimate center and radius at all, no interpretation.
        %Set up the return estimates
        center = cur_center.*chosen_mult*inv(rotate_matrix)+empirical_mean;
        %center = cur_center*inv(rotate_matrix)+empirical_mean;
        radius = cur_radius;
    %fprintf('radius and center')
    %radius
    %center
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Step 6: Project and rotate the sphere back into full space%
    
        augument_matrix = data_matrix.*0;
        augument_matrix(:,chosen_index) = X_rotated(:,chosen_index);
        augument_matrix = ProjectToSphere(augument_matrix,cur_center.*chosen_mult,cur_radius,2);

        % reduced matrix, please check

        % re-center at origin in the new coordinate frame
        reduced_matrix = augument_matrix - ones(size(data_matrix,1),1)*cur_center;
        % drop zero coordinates to reduce the dimension
        reduced_matrix = reduced_matrix(:,chosen_index);

        output_matrix = augument_matrix*inv(rotate_matrix); 
        output_matrix = output_matrix + empirical_mean;
        %disp(chosen_index)
        
        opt_ind = chosen_mult;
        %lossfun = @(x)ALG_LOSS(x,output_matrix,W_matrix);
        %[c1,r1] = EstimateSphere(output_matrix,lossfun)
    end
    %                                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
end
