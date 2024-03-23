function [best_vars, best_val,n_func_eval] = brute_force(func, bounds)
    % Generate all combinations of integer values within the bounds
    num_dims = size(bounds, 1);
    ranges = cell(1, num_dims);
    for i = 1:num_dims
        ranges{i} = bounds(i, 1):bounds(i, 2);
    end
    all_combinations = all_combinations_of(ranges);
    best_val = Inf;
    best_vars = [];
    count_c = 0;
    for i = 1:size(all_combinations, 1)
        vars = all_combinations(i, :);
        val = func(vars);
        count_c = count_c + 1;
        if val < best_val
            best_val = val;
            best_vars = vars;
        end
    end
    disp(['Brute-force: Number of evaluations = ', num2str(count_c)]);
    n_func_eval = count_c;
end