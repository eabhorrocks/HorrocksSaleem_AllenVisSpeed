function [padj_list, stat_list, df_list, stat_string_list] = extractPairwiseStats(varargin)
% extractPairwiseStats Extracts pairwise stats from the upper triangle of matrices.
% This function is flexible and accepts three input patterns:
%
% USAGE 1 (from multcompare-like stats struct):
%   [p, stat, df, str] = extractPairwiseStats(padj_matrix, stats_struct)
%
% USAGE 2 (from F-test matrices, e.g., coefTest):
%   [p, stat, df, str] = extractPairwiseStats(p_matrix, f_matrix, df1_matrix, df2_matrix)
%
% USAGE 3 (from Chi-squared matrices):
%   [p, stat, df, str] = extractPairwiseStats(p_matrix, x2_matrix, df_matrix)
%
% OUTPUTS:
%   padj_list        - A (N*(N-1)/2)x1 vector of p-values.
%   stat_list        - A (N*(N-1)/2)x1 vector of F-stats or X2-stats.
%   df_list          - A (N*(N-1)/2)xK vector of degrees of freedom (K=1 or 2).
%   stat_string_list - A (N*(N-1)/2)x1 cell array of formatted strings.

% --- Check which input pattern was used ---
if nargin == 2
    % --- CASE 1: (padj_matrix, stats_struct) ---
    padj = varargin{1};
    stats = varargin{2};
    
    if ~isnumeric(padj) || ~isstruct(stats)
        error('Usage 1: extractPairwiseStats(padj_matrix, stats_struct_matrix)');
    end
    
    num_items = size(padj, 1);
    num_pairs = (num_items * (num_items - 1)) / 2;
    
    % Pre-allocate outputs
    padj_list = nan(num_pairs, 1);
    stat_list = nan(num_pairs, 1);
    stat_string_list = cell(num_pairs, 1);
    
    % Try to get DF dimensions correctly
    try
        df_sample = stats(1, 2).df;
        df_cols = numel(df_sample);
        df_list = nan(num_pairs, df_cols);
    catch
        warning('Could not auto-detect df size. Assuming 2 columns.');
        df_list = nan(num_pairs, 2);
    end
    
    % Loop through the upper triangle
    k = 1; 
    for i = 1:(num_items - 1)
        for j = (i + 1):num_items
            padj_list(k) = padj(i, j);
            
            if isfield(stats(i,j), 'fstat') && ~isempty(stats(i,j).fstat)
                f_val = stats(i, j).fstat;
                df_vals = stats(i, j).df;
                
                stat_list(k) = f_val;
                df_list(k, :) = df_vals;
                
                % Create F-test string
                stat_string_list{k} = sprintf('F(%.0f, %.0f) = %.2f', df_vals(1), df_vals(2), f_val);
            else
                stat_string_list{k} = ''; % Add empty string if no stat
            end
            k = k + 1;
        end
    end

elseif nargin == 4
    % --- CASE 2: (p_matrix, f_matrix, df1_matrix, df2_matrix) ---
    padj  = varargin{1};
    Fstat = varargin{2};
    df1   = varargin{3};
    df2   = varargin{4};
    
    if ~isnumeric(padj) || ~isnumeric(Fstat) || ~isnumeric(df1) || ~isnumeric(df2)
        error('Usage 2: extractPairwiseStats(p_matrix, f_matrix, df1_matrix, df2_matrix)');
    end
    
    num_items = size(padj, 1);
    num_pairs = (num_items * (num_items - 1)) / 2;
    
    % Pre-allocate outputs
    padj_list = nan(num_pairs, 1);
    stat_list = nan(num_pairs, 1);
    df_list = nan(num_pairs, 2); 
    stat_string_list = cell(num_pairs, 1);

    % Loop through the upper triangle
    k = 1;
    for i = 1:(num_items - 1)
        for j = (i + 1):num_items
            f_val = Fstat(i, j);
            df1_val = df1(i, j);
            df2_val = df2(i, j);
            
            padj_list(k) = padj(i, j);
            stat_list(k) = f_val;
            df_list(k, :) = [df1_val, df2_val];
            
            % Create F-test string
            stat_string_list{k} = sprintf('F(%.0f, %.0f) = %.2f', df1_val, df2_val, f_val);
            
            k = k + 1;
        end
    end
    
elseif nargin == 3
    % --- CASE 3: (p_matrix, x2_matrix, df_matrix) ---
    padj   = varargin{1};
    X2stat = varargin{2};
    df     = varargin{3};
    
    if ~isnumeric(padj) || ~isnumeric(X2stat) || ~isnumeric(df)
        error('Usage 3: extractPairwiseStats(p_matrix, x2_matrix, df_matrix)');
    end
    
    num_items = size(padj, 1);
    num_pairs = (num_items * (num_items - 1)) / 2;
    
    % Pre-allocate outputs
    padj_list = nan(num_pairs, 1);
    stat_list = nan(num_pairs, 1); % (This will be X2)
    df_list = nan(num_pairs, 1);   % (X2 has 1 df)
    stat_string_list = cell(num_pairs, 1);

    % Loop through the upper triangle
    k = 1;
    for i = 1:(num_items - 1)
        for j = (i + 1):num_items
            p_val = padj(i, j);
            x2_val = X2stat(i, j);
            df_val = df(i, j);
            
            padj_list(k) = p_val;
            stat_list(k) = x2_val;
            df_list(k) = df_val;
            
            % Create Chi-squared string (using 'X2')
            stat_string_list{k} = sprintf('X2(%.0f) = %.2f', df_val, x2_val);
            
            k = k + 1;
        end
    end
    
else
    % --- ERROR: Invalid number of arguments ---
    error('Invalid number of arguments. Use 2, 3, or 4 inputs.');
end

end