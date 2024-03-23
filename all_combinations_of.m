function combs = all_combinations_of(ranges)
    % Generate all combinations from a cell array of ranges
    num_dims = numel(ranges);
    [grids{1:num_dims}] = ndgrid(ranges{:});
    combs = cell2mat(cellfun(@(c) c(:), grids, 'UniformOutput', false));
    combs = reshape(combs, [], num_dims);
end