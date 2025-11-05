function [F_string_list] = formatFStats(F_list, df1_list, df2_list)
% formatFStats Formats F-statistics and DFs into a cell array of strings.
%
% INPUTS:
%   F_list   - A vector (Nx1 or 1xN) of F-statistic values.
%   df1_list - A vector (same size as F_list) of numerator DFs.
%   df2_list - A vector (same size as F_list) of denominator DFs.
%
% OUTPUT:
%   F_string_list - An Nx1 cell array of formatted strings.

% Ensure all inputs are columns for consistent looping
F_list = F_list(:);
df1_list = df1_list(:);
df2_list = df2_list(:);

% Check for size mismatch
if numel(F_list) ~= numel(df1_list) || numel(F_list) ~= numel(df2_list)
    error('Input vectors (F, df1, df2) must all be the same size.');
end

num_items = numel(F_list);
F_string_list = cell(num_items, 1); % Pre-allocate cell array

% Loop and create the string for each item
for k = 1:num_items
    F_string_list{k} = sprintf('F(%.0f, %.0f) = %.2f', df1_list(k), df2_list(k), F_list(k));
end

end