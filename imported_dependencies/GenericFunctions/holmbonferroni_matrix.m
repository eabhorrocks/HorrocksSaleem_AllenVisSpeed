function [p_adj_matrix, h_matrix] = holmbonferroni_matrix(p_matrix, alpha)
%HOLMBONFERRONI_MATRIX Performs Holm-Bonferroni correction on a p-value matrix.
%
%   VERSION 2.0 - Corrected a bug that produced an all-NaN matrix.
%
%   DESCRIPTION:
%   This function is for a square, symmetric matrix of p-values (p_matrix)
%   from all-pairwise comparisons. It applies the Holm-Bonferroni correction
%   to the unique tests in the upper triangle.
%
%   INPUTS:
%   p_matrix    - An (n x n) symmetric matrix of p-values.
%   alpha       - (Optional) The significance level. Default is 0.05.
%
%   OUTPUTS:
%   p_adj_matrix - An (n x n) matrix of adjusted p-values. Diagonal is NaN.
%   h_matrix     - An (n x n) logical matrix indicating significance.

% --- Input Validation ---
if nargin < 1, error('Function requires a p-value matrix.'); end
if ~ismatrix(p_matrix) || size(p_matrix,1) ~= size(p_matrix,2)
    error('Input must be a square matrix.');
end
if nargin < 2 || isempty(alpha), alpha = 0.05; end

% --- Extract Unique P-values ---
n = size(p_matrix, 1);
upper_triangle_indices = find(triu(true(n), 1)); 
p_vector = p_matrix(upper_triangle_indices);
m = length(p_vector);

if m == 0
    p_adj_matrix = nan(n);
    h_matrix = false(n);
    return;
end

% --- Holm-Bonferroni Correction ---
[sorted_p, sort_idx] = sort(p_vector);
adj_p_sorted = zeros(size(sorted_p));
adj_p_sorted(1) = sorted_p(1) * m;
for i = 2:m
    adj_p_sorted(i) = max(adj_p_sorted(i-1), sorted_p(i) * (m - i + 1));
end
adj_p_sorted(adj_p_sorted > 1) = 1;
adj_p_vector = zeros(size(p_vector));
adj_p_vector(sort_idx) = adj_p_sorted;

% --- Reconstruct the Output Matrix ---
% *** FIX: Initialize with zeros, not NaNs ***
p_adj_matrix = zeros(n); 
p_adj_matrix(upper_triangle_indices) = adj_p_vector;

% Now, adding the transpose correctly fills the lower triangle
p_adj_matrix = p_adj_matrix + p_adj_matrix';

% *** FIX: Set the diagonal to NaN at the end ***
p_adj_matrix(logical(eye(n))) = NaN;

% Generate the logical significance matrix
h_matrix = p_adj_matrix < alpha;

end