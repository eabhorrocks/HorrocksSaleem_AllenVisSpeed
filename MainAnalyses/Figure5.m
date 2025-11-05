%% Figure 5

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

r2_thresh = 0.1;
r2p_thresh = 0.05;
rfp_thresh = 1;
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);

%% get tuning strength values for each area

for iarea = 1:8
    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));

    % stat
    allr2 = [areaUnits.r2dir_stat]';
    allr2p = [areaUnits.r2dirpval_stat]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),5,1);
    allTF = repmat([1:5]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,5,1);

    % remove instances where insufficient trials to calculate
    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allTF(idx)=[];
    allSigRF(idx)=[];

    ta(iarea).r2_stat = allr2;
    ta(iarea).allr2p_stat = allr2p;

    ta(iarea).areaVec_stat = allArea;
    ta(iarea).sessionVec_stat = allSession;
    ta(iarea).stateVec_stat = zeros(size(ta(iarea).sessionVec_stat));
    ta(iarea).tfVec_stat = allTF;
    ta(iarea).sigRFVec_stat = allSigRF;
    ta(iarea).nTuningCurves_stat = numel(allr2p);
    ta(iarea).nTuned_stat = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_stat = ta(iarea).nTuned_stat/ta(iarea).nTuningCurves_stat;
    ta(iarea).tunedFlag_stat = ta(iarea).r2_stat>r2_thresh & ta(iarea).allr2p_stat<r2p_thresh;


    % run
    allr2 = [areaUnits.r2dir_run]';
    allr2p = [areaUnits.r2dirpval_run]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),5,1);
    allTF = repmat([1:5]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,5,1);

    % remove instances where insufficient trials to calculate
    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allTF(idx)=[];
    allSigRF(idx)=[];

    ta(iarea).r2_run = allr2;
    ta(iarea).allr2p_run = allr2p;
    ta(iarea).areaVec_run = allArea;
    ta(iarea).sessionVec_run = allSession;
    ta(iarea).stateVec_run = ones(size(ta(iarea).sessionVec_run));
    ta(iarea).tfVec_run = allTF;
    ta(iarea).sigRFVec_run = allSigRF;
    ta(iarea).nTuningCurves_run = numel(allr2p);
    ta(iarea).nTuned_run = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_run = ta(iarea).nTuned_run/ta(iarea).nTuningCurves_run;
    ta(iarea).tunedFlag_run = ta(iarea).r2_run>r2_thresh & ta(iarea).allr2p_run<r2p_thresh;
    

    % combine stat and run
    ta(iarea).stateVec = categorical(cat(1,ta(iarea).stateVec_stat,ta(iarea).stateVec_run));
    ta(iarea).tunedVec = cat(1,ta(iarea).tunedFlag_stat,ta(iarea).tunedFlag_run);
    ta(iarea).r2vec = cat(1,ta(iarea).r2_stat, ta(iarea).r2_run);
    ta(iarea).sessionVec = categorical(cat(1,ta(iarea).sessionVec_stat,ta(iarea).sessionVec_run));
    ta(iarea).tfVec = categorical(cat(1,ta(iarea).tfVec_stat, ta(iarea).tfVec_run));
    ta(iarea).RFVec = categorical(cat(1,ta(iarea).sigRFVec_stat, ta(iarea).sigRFVec_run));
    ta(iarea).areaVec = categorical(cat(1,ta(iarea).areaVec_stat, ta(iarea).areaVec_run));

end


%% imagesc plot of p(tuned)
figure
imagesc([[ta.pTuned_stat]; [ta.pTuned_run]]), colormap(cmocean('tempo')), caxis([0.1 0.62])
ax=gca;ax.XTick=1:8; ax.XTickLabel=areas; colorbar
figure
imagesc([[ta.pTuned_run]-[ta.pTuned_stat]]), colormap(crameri('vik')), caxis([-0.22 0.22])
colorbar


%% glme on all units for p(tuned), to test effect of state

stateVec =cat(1,ta.stateVec);
tunedVec = cat(1,ta.tunedVec);
sessionVec = cat(1,ta.sessionVec);
tfVec = cat(1,ta.tfVec);
RFVec = cat(1,ta.RFVec);
areaVec = cat(1,ta.areaVec);

tbl = table(tunedVec, stateVec, sessionVec,tfVec,RFVec,areaVec,...
    'VariableNames', {'tunedFlag', 'state', 'session','tf','RF','area'});

f = 'tunedFlag ~ -1 + area:state'; % reduced models performs best
glme = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

for iarea = 1:8
    H=zeros(1,16);
    H(iarea*2-1)=-1;
    H(iarea*2)=1;
    H;
    [ta(iarea).pValue2, Fvals(iarea), df1(iarea), df2(iarea)] = coefTest(glme,H);
end

padj = holmbonferroni([ta.pValue2]);
[p_str, p_stars] = format_p_values(padj)

[F_string_list] = formatFStats(Fvals, df1, df2)
for i = 1:numel(F_string_list)
    fprintf('%s\n', F_string_list{i}); % Prints just the string
end

%% plot p(tuned) difference maps with significant tests for stationary and locomotion

stat_pTuned_diffArray = nan(8);
run_pTuned_diffArray = nan(8);
statP = nan(8);
runP = nan(8);

% stationary
for iarea1 = 1:8
    for iarea2 = 1:8
        
        % delta in p(tuned)
        stat_pTuned_diffArray(iarea1, iarea2) = ta(iarea1).pTuned_stat - ta(iarea2).pTuned_stat;
        run_pTuned_diffArray(iarea1, iarea2) = ta(iarea1).pTuned_run - ta(iarea2).pTuned_run;
        
        % F-tests
        H=zeros(1,16);
        H(iarea1*2-1)=-1;
        H(iarea2*2-1)=1;
        [statP(iarea1,iarea2), statF(iarea1,iarea2), statdf1(iarea1,iarea2), statdf2(iarea1,iarea2)] ...
            = coefTest(glme,H);

        H=zeros(1,16);
        H(iarea1*2)=-1;
        H(iarea2*2)=1;
        [runP(iarea1,iarea2), runF(iarea1,iarea2), rundf1(iarea1,iarea2), rundf2(iarea1,iarea2)] ...
            = coefTest(glme,H);

    end
end

% process p-values appropriately, including multiple comparisons correction
statPadj = holmbonferroni_matrix(statP);
[statP_str, statP_stars] = format_p_values(statPadj);

runPadj = holmbonferroni_matrix(runP);
[runP_str, runP_stars] = format_p_values(runPadj);


figure, hold on
imagesc(stat_pTuned_diffArray), axis ij, colormap(crameri('vik')), colorbar, caxis([-0.3 0.3])
for iarea1 = 1:8
    for iarea2 = 1:8
        text(iarea1,iarea2, statP_stars{iarea1,iarea2},'HorizontalAlignment', 'Center');
    end
end
title('stat p(tuned)')
hold on
for iarea = 1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
box on
xlim([0.5, 8.5]), ylim([0.5 8.5])
defaultAxesProperties(gca,false)

figure, hold on
imagesc(run_pTuned_diffArray), axis ij, axis ij, colormap(crameri('vik')), colorbar, caxis([-0.3 0.3])
for iarea1 = 1:8
    for iarea2 = 1:8
               text(iarea1,iarea2, runP_stars{iarea1,iarea2},'HorizontalAlignment', 'Center');

    end
end
title('run p(tuned)')
hold on
for iarea = 1:numel(areas)
    plot([iarea+0.5, iarea+0.5], [0 8.5],'k');
    plot([0 8.5], [iarea+0.5 iarea+0.5],'k')
end
ax=gca; ax.XTick = 1:8; ax.YTick = 1:8; ax.XTickLabels = areas; ax.YTickLabels = areas;
box on
xlim([0.5, 8.5]), ylim([0.5 8.5])
defaultAxesProperties(gca,false)


% format for stats sheet
% format data for google sheet stats
[padj_list, F_list, df_list,F_string_list] = extractPairwiseStats(statPadj, statF, statdf1, statdf2)

for i = 1:numel(F_string_list)
    fprintf('%s\n', F_string_list{i}); % Prints just the string
end

% format data for google sheet stats
[padj_list, F_list, df_list,F_string_list] = extractPairwiseStats(runPadj, runF, rundf1, rundf2)

for i = 1:numel(F_string_list)
    fprintf('%s\n', F_string_list{i}); % Prints just the string
end




%% Tuning strength survival functions and LME model

figure
for iarea = 1:8
    subplot(2,4,iarea), 
    hold on
    statvals = ta(iarea).r2_stat; statvals(statvals<0)=-.1;
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;

    [f, x] = ecdf(statvals);
    plot(x,1-f,'Color', areacols(iarea,:).*0.7);
    [f, x] = ecdf(runvals);
    plot(x,1-f,'Color', areacols(iarea,:));
    xlim([0 1])
    ylim([0 0.8])
    grid on
    plot([0.1 0.1], [0 1], 'k:')
    % xlim([])
    ax=gca; ax.YTick = 0:0.1:1; ax.XTick = 0:0.1:1;
    nvals(iarea,:) = [numel(statvals), numel(runvals)];
    defaultAxesProperties(gca,true)
end

stateVec =cat(1,ta.stateVec);
r2Vec = cat(1,ta.r2vec); r2Vec(r2Vec<0)=0;
sessionVec = cat(1,ta.sessionVec);
tfVec = cat(1,ta.tfVec);
RFVec = cat(1,ta.RFVec);
areaVec = cat(1,ta.areaVec);

tbl = table(r2Vec, stateVec, sessionVec,tfVec,RFVec,areaVec,...
    'VariableNames', {'tuningStrength', 'state', 'session','tf','RF','area'});

f = 'tuningStrength ~ -1 + area:state + (1|session) + (1|tf) + (1|RF)';
%f = 'tunedFlag ~ -1 + area:state + (1|dir) + (1|RF)';

lme = fitlme(tbl,f,'DummyVarCoding','full');

H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1*2-1) = -1; H_temp(iarea1*2) = 1;
    [pvals(iarea1), Fvals(iarea1), df1(iarea1), df2(iarea1)] = coefTest(lme,H_temp);
end

padj = holmbonferroni(pvals);
[p_str, p_stars] = format_p_values(padj)

[F_string_list] = formatFStats(Fvals, df1, df2)
for i = 1:numel(F_string_list)
    fprintf('%s\n', F_string_list{i}); % Prints just the string
end

%% example tuning curves

IDs = [950935552, 951769160, 951855640];
itf = [1, 3 1];


for iunit = 1:numel(IDs)
    subplot(1,3,iunit), hold on

    idx = find([goodUnits.ID]==IDs(iunit));

     hold on
    for idir = 1:8

        plot(repelem(idir,1,numel(goodUnits(idx).spikeCounts_run{idir, itf(iunit)})),...
            goodUnits(idx).spikeCounts_run{idir,itf(iunit)}, 'k.', 'MarkerSize',10)
    end
    
    plot(1:8, goodUnits(idx).tuning_run(:,itf(iunit)), 'k')
    title(goodUnits(idx).r2dir_run(itf(iunit)))


    % plot(1:0.01:7, feval(gaussFun, areaUnits(idx(iunit)).gaussParams_run(idirs(iunit),:),  1:0.01:7), 'r-')
    % % ax = gca; ax.XTick = 1:7;
    % defaultAxesProperties(gca, true)

    % subplot(122)
    % plot(1:7, areaUnits(idx(iunit)).SSI_run(idirs(iunit),:), 'k-')
    % ylim([0.5 1.9]); 
    % % ax= gca; ax.YTick = 0.5:0.2:1.9; ax.XTick = 1:7;
    defaultAxesProperties(gca, true)

end
