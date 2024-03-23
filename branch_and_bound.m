function [best_vars, best_val] = branch_and_bound(func, bounds, tolerance)
    queue = {bounds};
    best_val = inf;
    best_vars = [];
    n_func_eval = 0;

    while ~isempty(queue)
        current_bounds = queue{1};
        queue(1) = [];
function [best_vars, best_val,n_func_eval] = branch_and_bound(func, bounds, tolerance)
    queue = {bounds};
    best_val = inf;
    best_vars = [];
    n_func_eval = 0;

    while ~isempty(queue)
        current_bounds = queue{1};
        queue(1) = [];

        % Calculate the midpoint of each dimension in current_bounds
        mid_vars = zeros(1, size(current_bounds, 1));
        for i = 1:size(current_bounds, 1)
            mid_vars(i) = floor((current_bounds(i, 1) + current_bounds(i, 2)) / 2);
        end

        val = func(mid_vars);
        n_func_eval = n_func_eval + 1;
        if val < best_val
            best_val = val;
            best_vars = mid_vars;
        end

        % Split bounds if larger than tolerance
        for i = 1:size(current_bounds, 1)
            if current_bounds(i, 2) - current_bounds(i, 1) > tolerance
                mid_point = floor((current_bounds(i, 1) + current_bounds(i, 2)) / 2);
                left_bounds = current_bounds;
                right_bounds = current_bounds;
                left_bounds(i, 2) = mid_point;
                right_bounds(i, 1) = mid_point + 1;
                queue{end + 1} = left_bounds;
                queue{end + 1} = right_bounds;
                break;
            end
        end
    end
    disp(['Branch-and-Bound: Number of evaluations = ', num2str(n_func_eval)]);
end
 
        % Evaluate at the midpoint
        mid_vars = arrayfun(@(b) (b(1) + b(2)) // 2, current_bounds);
        val = func(mid_vars);
        n_func_eval = n_func_eval + 1;
        if val < best_val
            best_val = val;
            best_vars = mid_vars;
        end

        % Split bounds if larger than tolerance
        for i = 1:size(current_bounds, 1)
            if current_bounds(i, 2) - current_bounds(i, 1) > tolerance
                mid_point = (current_bounds(i, 1) + current_bounds(i, 2)) // 2;
                left_bounds = current_bounds;
                right_bounds = current_bounds;
                left_bounds(i, 2) = mid_point;
                right_bounds(i, 1) = mid_point + 1;
                queue{end + 1} = left_bounds;
                queue{end + 1} = right_bounds;
                break;
            end
        end
    end
    disp(['Branch-and-Bound: Number of evaluations = ', num2str(n_func_eval)]);
end