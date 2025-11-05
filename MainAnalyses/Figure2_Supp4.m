%% Visual speed tuning class by brain areacharacterisation
% Figure 2
% addpath 'D:\Code\GenericFunctions'
load('goodUnits.mat')

%%
areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
r2_thresh = 0.1;
r2p_thresh=0.05;
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);

%original: 1 = lowpass, 2 = highpass, 3 = bandpass, 4 = inverted
%plot order: 1 = lowpass, 2 = bandpass, 3 = inverted, 4 = highpass
plotOrder = [1 3 4 2];
%% get distributions of preferred speed and shape by visual area

for iarea = 1:numel(areas)

    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));

    % tuned stat units
    idx = find([areaUnits.r2_stat]>=r2_thresh & [areaUnits.r2pval_stat]<=r2p_thresh);

    allPrefSpeeds = [areaUnits.prefSpeed_stat];
    allGaussChar = [areaUnits.gaussChar_stat];
    allGaussR2 = [areaUnits.gaussR2_stat];
    allr2 = [areaUnits.r2_stat];
    allr2p = [areaUnits.r2pval_stat];
    allSession = repelem([areaUnits.sessionID],1,4);
    allDir = repmat(1:4,1,numel(areaUnits));

    statta(iarea).prefSpeeds = allPrefSpeeds(idx);
    statta(iarea).gaussChar = allGaussChar(idx);
    statta(iarea).gaussR2 = allGaussR2(idx);
    statta(iarea).r2 = allr2(idx);
    statta(iarea).r2p = allr2p(idx);
    statta(iarea).session = allSession(idx);
    statta(iarea).dir = allDir(idx);

    % tuned stat units
    idx = find([areaUnits.r2_run]>=r2_thresh & [areaUnits.r2pval_run]<=r2p_thresh);

    allPrefSpeeds = [areaUnits.prefSpeed_run];
    allGaussChar = [areaUnits.gaussChar_run];
    allGaussR2 = [areaUnits.gaussR2_run];
    allr2 = [areaUnits.r2_run];
    allr2p = [areaUnits.r2pval_run];

    runta(iarea).prefSpeeds = allPrefSpeeds(idx);
    runta(iarea).gaussChar = allGaussChar(idx);
    runta(iarea).gaussR2 = allGaussR2(idx);
    runta(iarea).r2 = allr2(idx);
    runta(iarea).r2p = allr2p(idx);
    runta(iarea).session = allSession(idx);
    runta(iarea).dir = allDir(idx);
end

allPrefSpeeds = cat(2,statta.prefSpeeds,runta.prefSpeeds)';
allGaussChar = categorical(cat(2,statta.gaussChar,runta.gaussChar)');
allAreaVec = categorical(cat(2,repelem(1:8,1,cellfun(@numel, {statta.prefSpeeds})),repelem(1:8,1,cellfun(@numel, {runta.prefSpeeds}))))';
allDir = categorical(cat(2,statta.dir,runta.dir)');
allSession = categorical(cat(2,statta.session,runta.session)');

%% tuning class dists by area

gcVals= nan(4,8);

for iarea = 1:8
    for igc = 1:4
        gcVals(igc,iarea)=numel(find([[statta(iarea).gaussChar], [runta(iarea).gaussChar]]==igc));
    end
end

gcVals_norm = gcVals./sum(gcVals,1);
gcVals_norm = gcVals_norm(plotOrder,:);

figure
bar(gcVals_norm', 'stacked')
defaultAxesProperties(gca, true)
ax=gca; ax.XTick = 1:8; ax.XTickLabel =areas;
ylim([0 1]), legend({'lowpass','banpdass', 'bandreject','highpass'})

% X2 test for tuning class and brain area
GCvec = double(allGaussChar);
areavec = double(allAreaVec);
countTable = nan(8,4);
for iarea = 1:8
    for igc = 1:4
        countTable(iarea,igc) = sum([GCvec]==igc & [areavec]==iarea);
    end
end

[X2stat, p, df] = chi2test_countTable(countTable)

% pairwise comparisons between areas
pArray = nan(8);
for iarea1 = 1:8
    for iarea2 = 1:8
        [X2array(iarea1,iarea2), pArray(iarea1,iarea2), df(iarea1,iarea2)] = chi2test_countTable(countTable([iarea1 iarea2],:));
    end
end
padj = holmbonferroni_matrix(pArray);
[p_str, p_stars] = format_p_values(padj);
padj(1:size(padj,1)+1:end) = 1;

figure
tickvals = [0 0.001 0.01 0.05 1];
pvals_disc = discretize(padj,tickvals,'IncludedEdge','right');
colmap = colormap(slanCM('gray',6));
colmap=flipud(colmap);
colmap = colmap(1:4,:);
colmap = [colmap; 0 0 0];
imagesc(pvals_disc), colormap(colmap); 
c = colorbar;
c.Ticks =linspace(1,4,6); c.TickLabels = {"","0.0001","0.001", "0.01", "0.05", ">0.05"};

ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
defaultAxesProperties(gca,false)
hold on

for iarea = 1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end

[p_list, X2_list, df_list, str_list] = extractPairwiseStats(padj, X2array, df);
for i = 1:numel(str_list)
    fprintf('%s\n', str_list{i}); % Prints just the string
end


%% Preferred speed distributions

figure
for iarea = 1:8
    subplot(2,4,iarea)
    ta(iarea).prefSpeeds = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    h = histcounts(ta(iarea).prefSpeeds);
    bar(1:7,h./sum(h), 'FaceColor', areacols(iarea,:))
    hold on
    ylim([0 0.5]); ax = gca; ax.XTick = 1:7; ax.YTick = 0:0.1:0.5;
    ta(iarea).rawmean = mean(ta(iarea).prefSpeeds);
    ta(iarea).meanPref = 2^(3+mean([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]]));
    ta(iarea).varPref = var([[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]]);
    meanval = mean(ta(iarea).prefSpeeds);
    plot([meanval-ta(iarea).varPref/2, meanval+ta(iarea).varPref/2],[0.45, 0.45],'Color',areacols(iarea,:))
    plot(mean(ta(iarea).prefSpeeds), 0.45,'o','MarkerFaceColor',areacols(iarea,:),'MarkerEdgeColor','w')
    title(areas{iarea})
    nvals(iarea) = sum(h);
    defaultAxesProperties(gca, false)
end

%% Centre of mass of preferred speeds 

figure, hold on
for iarea = 1:8
    vals = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    mean_val = mean(vals);
    sem_val = sem(vals(:)).*1.96;
    plot(iarea,mean_val, 'o', 'MarkerFaceColor', areacols(iarea,:), 'MarkerEdgeColor', 'w')
    plot([iarea, iarea], [mean_val-sem_val, mean_val+sem_val], 'Color', areacols(iarea,:));

end
ylim([3 5])
ax=gca; ax.XTick = 1:8; ax.XTickLabel = areas;
ax.YTick = [3 4 5]; ax.YTickLabel = [32 64 128]; ylabel('Visual speed (deg/s)');
defaultAxesProperties(gca, true)

%linear mixed effects model for statistical analysis of preferred speeds
tbl = table(allPrefSpeeds, allGaussChar, allAreaVec,allDir, allSession,...
    'VariableNames', {'prefSpeeds', 'gaussChar', 'areas', 'dir', 'session'});

f = 'prefSpeeds ~ -1 + areas + (1|session) + (1|dir)';
lme = fitlme(tbl,f,'dummyVarCoding','full')

% pairwise F-tests
for iarea1 = 1:8
    for iarea2 = 1:8
        H=zeros(1,8);
        H(iarea1)=1;H(iarea2)=-1;
        [pvals(iarea1,iarea2), F(iarea1,iarea2), df1(iarea1,iarea2), df2(iarea1,iarea2)]...
            = coefTest(lme,H);
    end
end

% mult comparisons correction
padj = holmbonferroni_matrix(pvals);
[p_str, p_stars] = format_p_values(padj);

% calculate difference in means and plot using imagesc
for iarea1 = 1:8
    for iarea2 = 1:8
        valDiff(iarea1,iarea2) = log2(ta(iarea1).meanPref)-log2(ta(iarea2).meanPref);
    end
end

figure, hold on
imagesc(valDiff), axis ij
for iarea1 = 1:8
    for iarea2 = 1:8
        text(iarea1,iarea2, p_stars{iarea1,iarea2},'HorizontalAlignment', 'Center');
    end
end
colormap(crameri('vik')), colorbar

% plot grids
hold on
for iarea = 1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
box on
xlim([0.5, 8.5]), ylim([0.5 8.5])
defaultAxesProperties(gca,false)


% format data for google sheet stats
[padj_list, F_list, df_list,F_string_list] = extractPairwiseStats(padj, F, df1, df2)

for i = 1:numel(F_string_list)
    fprintf('%s\n', F_string_list{i}); % Prints just the string
end





%% Variance of preferred speed distributions

figure, hold on
for iarea = 1:8
    vals = [[statta(iarea).prefSpeeds],[runta(iarea).prefSpeeds]];
    var_val = var(vals);
    [h,p,ci,stats] = vartest(vals,0, 'Alpha', 0.05);
    plot(iarea,var_val, 'o', 'MarkerFaceColor', areacols(iarea,:), 'MarkerEdgeColor', 'w')
    plot([iarea, iarea], [ci(1), ci(2)], 'Color', areacols(iarea,:));

end
ylim([1 4]); ax=gca; ax.XTick=1:8; ax.XTickLabel = areas; ylabel('Variance (log_2 deg/s)')
defaultAxesProperties(gca, true)

% statistical test of variance of speed distributions
for iarea = 1:8
    ta(iarea).prefSpeeds = cat(2,statta(iarea).prefSpeeds, runta(iarea).prefSpeeds);
end

group = repelem(1:8, [cellfun(@numel, {ta.prefSpeeds})]);
vals = [ta.prefSpeeds];

[p, stats] = vartestn(vals(:),group(:), 'TestType', 'LeveneQuadratic','display','off')

for iarea1 = 1:8
    for iarea2 = 1:8
        if iarea1==iarea2
            pArray(iarea1,iarea2) = nan;
        else
            group = repelem([iarea1, iarea2], [cellfun(@numel, {ta([iarea1, iarea2]).prefSpeeds})]);
            vals = [ta([iarea1, iarea2]).prefSpeeds];
            [pArray(iarea1,iarea2), stats(iarea1,iarea2)] = vartestn(vals(:),group(:), 'TestType', 'LeveneQuadratic','display', 'off');
        end
    end
end

for iarea1 = 1:8
    for iarea2 = 1:8
        valDiff(iarea1,iarea2) = ta(iarea1).varPref-ta(iarea2).varPref;
    end
end

padj = holmbonferroni_matrix(pArray);
[p_str, p_stars] = format_p_values(padj);
padj(1:size(padj,1)+1:end) = 1;

padj(padj>1)=1;
figure, hold on
imagesc(valDiff), axis ij
for iarea1 = 1:8
    for iarea2 = 1:8
        text(iarea1,iarea2, p_stars{iarea1,iarea2},'HorizontalAlignment', 'Center');
    end
end

colormap(crameri('vik')), colorbar
xlim([0.5 8.5]), ylim([0.5 8.5])
for iarea = 1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;

defaultAxesProperties(gca, false)

[padj_list, F_list, df_list, F_string_list] = extractPairwiseStats(p_adj, stats)

% get pairwise test results in order for google sheet
num_pairs = (numel(areas) * (numel(areas)-1))/2;
pairwise_list_padj = nan(num_pairs, 1);
pairwise_list_F = nan(num_pairs, 1);


% Loop through the upper triangle row-by-row
k = 1; % This is our list index
for iarea1 = 1:numel(areas)-1
    for iarea2 = (iarea1 + 1):numel(areas)
        % This loop gets M(1,2), M(1,3)...M(1,8), then
        % M(2,3), M(2,4)...M(2,8), etc.
        pairwise_list_padj(k) = padj(iarea1, iarea2);
        pairwise_list_F(k) = stats(iarea1, iarea2).fstat;
        pairwise_list_df(k,:) = stats(iarea1,iarea2).df;
        k = k + 1;
    end
end

