function [center,radius] = EstimateSphere(data_matrix,lossfun)
  %Function: Estimate the center and radius of the same dimensional sphere
  %as the given data matrix by minimizing the given lossfun
    data_dim = size(data_matrix);
    fun = lossfun;
    %nonlcon = confungrad;
    x0 = ones(1,data_dim(2)+1);
    %A = [];
    %b = [];
    %Aeq = [];
    %beq = [];
    %lb = [];
    %ub = [];
    options = optimoptions('fminunc',...
                            'Algorithm','quasi-newton',...
                            'Display','final',...
                            'MaxIterations',10000,...
                            'MaxFunctionEvaluations',10000,...
                            'OptimalityTolerance', 1e-12, ...
                            'StepTolerance', 1e-12, ...
                            'CheckGradients',true);
    %[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    [xval,fval] = fminunc(fun,x0,options);
    center = xval( 1:(length(xval)-1) );
    radius = xval( length(xval) );
end