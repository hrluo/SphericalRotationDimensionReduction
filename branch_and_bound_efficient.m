function [best_vars, best_val] = branch_and_bound_efficient(func, bounds, tolerance)
    queue = {bounds};
    best_val = inf;
    best_vars = [];
    n_func_eval = 0;
    evaluated_points = {}; % Cell array to store evaluated points

    while ~isempty(queue)
        current_bounds = queue{1};
        queue(1) = [];
function [best_vars, best_val,n_func_eval] = branch_and_bound_efficient(func, bounds, tolerance)
    queue = {bounds};
    best_val = inf;
    best_vars = [];
    n_func_eval = 0;
    evaluated_points = {}; % Cell array to store evaluated points as strings

    while ~isempty(queue)
        current_bounds = queue{1};
        queue(1) = [];

        % Calculate the midpoint of each dimension in current_bounds
        mid_vars = zeros(1, size(current_bounds, 1));
        for i = 1:size(current_bounds, 1)
            mid_vars(i) = floor((current_bounds(i, 1) + current_bounds(i, 2)) / 2);
        end

        mid_vars_str = mat2str(mid_vars); % Convert to string for comparison
        if ~ismember(mid_vars_str, evaluated_points)
            val = func(mid_vars);
            n_func_eval = n_func_eval + 1;
            evaluated_points{end + 1} = mid_vars_str; % Store as string
            if val < best_val
                best_val = val;
                best_vars = mid_vars;
            end
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
    disp(['Branch-and-Bound (Efficient): Number of evaluations = ', num2str(n_func_eval)]);
end

 
        % Evaluate at the midpoint
        mid_vars = arrayfun(@(b) (b(1) + b(2)) // 2, current_bounds);
        if ~ismember(mid_vars, evaluated_points, 'rows')
            val = func(mid_vars);
            n_func_eval = n_func_eval + 1;
            evaluated_points{end + 1} = mid_vars;
            if val < best_val
                best_val = val;
                best_vars = mid_vars;
            end
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
    disp(['Branch-and-Bound (Efficient): Number of evaluations = ', num2str(n_func_eval)]);
end