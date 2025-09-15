function [p_str, stars] = format_p_values(pvals)
%FORMAT_P_VALUES Converts numeric p-values to manuscript-ready strings and asterisks.
%
%   VERSION 2.0
%
%   SYNTAX:
%   [p_str, stars] = format_p_values(pvals)
%
%   DESCRIPTION:
%   Converts a numeric vector of p-values (pvals) into a cell array of
%   formatted strings (p_str) and a cell array of significance
%   asterisks (stars).
%
%   INPUT:
%   pvals   - A numeric vector of raw p-values.
%
%   OUTPUTS:
%   p_str   - A cell array of strings with formatted p-values (e.g., "p = 0.023").
%   stars   - A cell array of strings with significance stars:
%             '****' for p < 0.0001
%             '***'  for p < 0.001
%             '**'   for p < 0.01
%             '*'    for p < 0.05
%             ''     for p >= 0.05 (not significant)
%
%   EXAMPLE:
%   raw_p = [0.25, 0.045, 0.0023, 0.0005, 0.000012];
%   [p_string, p_stars] = format_p_values(raw_p);
%   T = table(raw_p(:), p_string(:), p_stars(:), ...
%       'VariableNames', {'RawP', 'Formatted', 'Stars'});
%   disp(T);

% --- Input Validation ---
if nargin < 1 || ~isnumeric(pvals)
    error('Function requires a numeric vector of p-values as input.');
end

% Initialize output cell arrays
p_str = cell(size(pvals));
stars = cell(size(pvals));

% --- Loop through each p-value and apply formatting ---
for i = 1:numel(pvals)
    p = pvals(i); % Get the current p-value
    
    % --- Part 1: Format the p-value string ---
    if p >= 0.001
        p_str{i} = sprintf('p = %.3f', p);
    elseif p >= 1e-4
        p_str{i} = 'p < 0.001';
    else
        exponent = ceil(-log10(p));
        p_str{i} = sprintf('p < 10^-%d', exponent);
    end
    
    % --- Part 2: Determine the significance stars ---
    if p < 0.0001
        stars{i} = '****';
    elseif p < 0.001
        stars{i} = '***';
    elseif p < 0.01
        stars{i} = '**';
    elseif p < 0.05
        stars{i} = '*';
    else
        stars{i} = 'n.s.'; % Assign an empty string for non-significant
    end
end

end