%% Supp Figure 5: Compare tuning strength with reliability and mutual informatin
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);
areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

%% calculate reliability
for iunit =1:numel(goodUnits)

    % calculate reliability from Christensen & Pillow 2022
    for idir = 1:4
        goodUnits(iunit).reli_stat(1,idir) = calcReliability(goodUnits(iunit).spikeCounts_stat(idir,:));
        goodUnits(iunit).reli_run(1,idir) = calcReliability(goodUnits(iunit).spikeCounts_run(idir,:));
    end
end

%%

for iarea = 1:8

    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));
    allr2_stat = cat(1,areaUnits.r2_stat);
    allMI_stat = cat(1,areaUnits.MI_stat);
    allreli_stat = cat(1,areaUnits.reli_stat);

    allr2_run = cat(1,areaUnits.r2_run);
    allMI_run = cat(1,areaUnits.MI_run);
    allreli_run = cat(1,areaUnits.reli_run);

    % Concatenate each pair vertically
    allr2    = [allr2_stat(:);    allr2_run(:)];
    allMI    = [allMI_stat(:);    allMI_run(:)];
    allreli  = [allreli_stat(:);  allreli_run(:)];

    % State variable: 0 for stat, 1 for run
    allstate = [zeros(size(allr2_stat(:))); ones(size(allr2_run(:)))];

    % Make a table
    areaTable = table(allr2, allMI, allreli, allstate);

    % Reliability vs R2
    [axScatter, axTop, axRight] = scatterWithMarginals(areaTable, 'allreli', 'allr2', 'allstate', [areacols(iarea,:)*0.7; areacols(iarea,:)], 0.6);
    axScatter.XLim = [0 1]; axScatter.YLim = [-0.5 1];  axScatter.YTick = -0.4:0.2:1;
    axScatter.XLabel.String = 'Reliability';  axScatter.YLabel.String = 'Tuning Strength (R^2)'; 
    defaultAxesProperties(axScatter,false)
    % get correlation and add to axes
    [rho,pval] = corr(areaTable.allreli,areaTable.allr2, Rows="complete");
    % Format rho
    rho_str = sprintf('%.2g', rho);
    % Format p-value in scientific notation if very small
    if pval < 0.001, p_str = 'p<0.001'; else, p_str = sprintf('p=%.3f', pval); end
    labelStr = ['r=' rho_str ', ' p_str];
    x = axScatter.XLim(1) + 0.02*range(axScatter.XLim);  % slightly right from left
    y = axScatter.YLim(2) - 0.02*range(axScatter.YLim);  % slightly below top
    text(axScatter, x, y, labelStr, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight','normal');
    print(gcf,[areas{iarea},'_','reli_vs_r2'],'-dsvg', '-vector');
    
    % MI vs R2
    [axScatter, axTop, axRight] = scatterWithMarginals(areaTable, 'allMI', 'allr2', 'allstate', [areacols(iarea,:)*0.7; areacols(iarea,:)], 0.6);
    axScatter.YLim = [-0.5 1];  axScatter.YTick = -0.4:0.2:1;
    axScatter.XLabel.String = 'Mutual Information (bits)';  axScatter.YLabel.String = 'Tuning Strength (R^2)'; 
    defaultAxesProperties(axScatter,false)
    % get correlation and add to axes
    [rho,pval] = corr(areaTable.allMI,areaTable.allr2, Rows="complete");
    % Format rho
    rho_str = sprintf('%.2g', rho);
    % Format p-value in scientific notation if very small
    if pval < 0.001, p_str = 'p<0.001'; else, p_str = sprintf('p=%.3f', pval); end
    labelStr = ['r=' rho_str ', ' p_str];
    x = axScatter.XLim(1) + 0.02*range(axScatter.XLim);  % slightly right from left
    y = axScatter.YLim(2) - 0.02*range(axScatter.YLim);  % slightly below top
    text(axScatter, x, y, labelStr, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight','normal');
    print(gcf,[areas{iarea},'_','MI_vs_r2'],'-dsvg', '-vector');
end

