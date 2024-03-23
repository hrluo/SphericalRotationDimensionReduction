function [best_vars, best_val, n_evals] = hl_minimize(func, bounds, method, tolerance)
    if nargin < 4
        tolerance = 1; % Default tolerance
    end
    if nargin < 3
        method = 'branch_and_bound'; % Default method
    end

    if strcmp(method, 'brute_force')
        [best_vars, best_val, n_evals] = brute_force(func, bounds);
    elseif strcmp(method, 'branch_and_bound')
        [best_vars, best_val, n_evals] = branch_and_bound(func, bounds, tolerance);
    elseif strcmp(method, 'branch_and_bound_efficient')
        [best_vars, best_val, n_evals] = branch_and_bound_efficient(func, bounds, tolerance);
    else
        error('Invalid method specified. Choose "brute_force" or "branch_and_bound".');
    end
end