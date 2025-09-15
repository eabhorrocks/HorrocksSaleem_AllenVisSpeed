function [axScatter, axTop, axRight] = scatterWithMarginals(T, xVar, yVar, groupVar, groupColors, scatterFrac)
% SCATTERWITHMARGINALS Scatter plot with smoothed marginal histograms
%
% Returns axes handles: [axScatter, axTop, axRight]
%
% T: table containing your data
% xVar, yVar: names of variables to plot (strings)
% groupVar: name of grouping variable (optional)
% groupColors: nGroups x 3 RGB array for colors (optional)
% scatterFrac: fraction of figure width/height for main scatter panel (default 0.7)

if nargin < 4 || isempty(groupVar)
    T.Group = ones(height(T),1);
    groupVar = 'Group';
end

groups = unique(T.(groupVar));
nGroups = numel(groups);

if nargin < 5 || isempty(groupColors)
    groupColors = lines(nGroups); 
end

if nargin < 6 || isempty(scatterFrac)
    scatterFrac = 0.7; 
end

figure;

margin = 0.05; % spacing between axes

% --- Scatter panel position ---
scatterLeft = 0.1;
scatterBottom = 0.1;
scatterWidth = scatterFrac;
scatterHeight = scatterFrac;

% --- Main scatter plot ---
axScatter = axes('Position',[scatterLeft, scatterBottom, scatterWidth, scatterHeight]); hold on;
for i = 1:nGroups
    idx = T.(groupVar) == groups(i);
    scatter(axScatter, T.(xVar)(idx), T.(yVar)(idx), 36, groupColors(i,:), 'filled', 'MarkerFaceAlpha',0.6);
end
xlabel(axScatter, xVar); ylabel(axScatter, yVar);
axScatter.Box = 'off';      % remove box
axScatter.XGrid = 'off';    % remove x grid
axScatter.YGrid = 'off';    % remove y grid

% --- Top histogram ---
topHeight = 1 - scatterBottom - scatterHeight - margin;
axTop = axes('Position',[scatterLeft, scatterBottom + scatterHeight + margin, scatterWidth, topHeight]); hold on;
for i = 1:nGroups
    idx = T.(groupVar) == groups(i);
    [f,xi] = ksdensity(T.(xVar)(idx));
    plot(axTop, xi, f, 'LineWidth', 1.5, 'Color', groupColors(i,:));
end
axTop.XAxis.Visible = 'off';
axTop.YAxis.Visible = 'off';
axTop.Box = 'off';
axTop.Color = 'none';
ylim(axTop, 'auto'); % auto scale y-axis
axTop.XLim = axScatter.XLim;

% --- Right histogram ---
rightWidth = 1 - scatterLeft - scatterWidth - margin;
axRight = axes('Position',[scatterLeft + scatterWidth + margin, scatterBottom, rightWidth, scatterHeight]); hold on;
for i = 1:nGroups
    idx = T.(groupVar) == groups(i);
    [f,xi] = ksdensity(T.(yVar)(idx));
    plot(axRight, f, xi, 'LineWidth', 1.5, 'Color', groupColors(i,:));
end
axRight.XAxis.Visible = 'off';
axRight.YAxis.Visible = 'off';
axRight.Box = 'off';
axRight.Color = 'none';
xlim(axRight, 'auto'); % auto scale x-axis
axRight.YLim = axScatter.YLim;

% --- Link axes ---
linkaxes([axScatter, axTop],'x');
linkaxes([axScatter, axRight],'y');

end
