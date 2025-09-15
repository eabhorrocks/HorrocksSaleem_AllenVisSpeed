function [p_adj, h] = holmbonferroni(pvals, alpha)
%HOLMBONFERRONI Perform the Holm-Bonferroni correction for multiple comparisons.
%
%   SYNTAX:
%   [p_adj, h] = holmbonferroni(pvals)
%   [p_adj, h] = holmbonferroni(pvals, alpha)
%
%   DESCRIPTION:
%   This function computes adjusted p-values from a vector of raw p-values
%   (pvals) using the Holm-Bonferroni sequential method, which provides strong
%   control of the family-wise error rate (FWER) and is more powerful
%   than the standard Bonferroni correction.
%
%   INPUTS:
%   pvals       - A vector of raw p-values from your multiple tests.
%   alpha       - (Optional) The significance level to use. Default is 0.05.
%
%   OUTPUTS:
%   p_adj       - A vector of the same size as pvals containing the
%                 Holm-Bonferroni adjusted p-values.
%   h           - A logical vector of the same size as pvals. An element
%                 is true (1) if the corresponding adjusted p-value is
%                 less than alpha, indicating a significant result.
%
%   EXAMPLE:
%   % Suppose you have 8 p-values from an experiment
%   raw_p = [0.002, 0.045, 0.011, 0.230, 0.049, 0.001, 0.150, 0.021];
%   [adjusted_p, significant_results] = holmbonferroni(raw_p, 0.05);
%   % Display results in a table
%   T = table(raw_p', adjusted_p', significant_results', ...
%       'VariableNames', {'RawP', 'AdjustedP', 'IsSignificant'});
%   disp(T);

% --- Input Validation ---
if nargin < 1
    error('Function requires at least one input: a vector of p-values.');
end
if ~isvector(pvals) || ~isnumeric(pvals)
    error('Input pvals must be a numeric vector.');
end
if nargin < 2 || isempty(alpha)
    alpha = 0.05; % Default significance level
end

% --- Holm-Bonferroni Calculation ---
m = length(pvals); % Number of tests

% Sort p-values in ascending order, but keep track of their original indices
[sorted_p, sort_idx] = sort(pvals);

% Calculate the adjusted p-values using the sequential formula
adj_p_sorted = zeros(1, m);
adj_p_sorted(1) = sorted_p(1) * m;
for i = 2:m
    % The i-th adjusted p-value is the maximum of the previous one
    % and the current p-value multiplied by the number of remaining tests.
    adj_p_sorted(i) = max(adj_p_sorted(i-1), sorted_p(i) * (m - i + 1));
end

% Ensure adjusted p-values do not exceed 1
adj_p_sorted(adj_p_sorted > 1) = 1;

% Reorder the adjusted p-values to match the original input order
p_adj = zeros(1, m);
p_adj(sort_idx) = adj_p_sorted;

% Determine which results are significant
if nargout > 1
    h = p_adj < alpha;
end

end