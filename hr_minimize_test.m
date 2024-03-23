% Test case for hl_minimize and brute_force
example_func = @(vars) (vars(1) - 1)^2 + abs(vars(2) - 2) + vars(1)*vars(3) + (vars(4)^3 - vars(4) + 1);
bounds = [0, 9; -5, -1; 0, 9; 0, 19];

[best_vars, best_val, n_evals] = hl_minimize(example_func, bounds, 'brute_force');
disp(['Best variables: ', mat2str(best_vars)]);
disp(['Best value: ', num2str(best_val)]);

[best_vars, best_val, n_evals] = hl_minimize(example_func, bounds, 'branch_and_bound',1);
disp(['Best variables: ', mat2str(best_vars)]);
disp(['Best value: ', num2str(best_val)]);

[best_vars, best_val, n_evals] = hl_minimize(example_func, bounds, 'branch_and_bound_efficient',1);
disp(['Best variables: ', mat2str(best_vars)]);
disp(['Best value: ', num2str(best_val)]);