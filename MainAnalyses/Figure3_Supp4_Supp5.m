%% Comparison of tuning strength between states in each area
% Figure 3

addpath 'D:\Code\GenericFunctions'
load('goodUnits.mat')

%%
r2_thresh = 0.1; r2p_thresh = 0.05;
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);
areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

%% P(Tuned)

% find fraction tuned for each area in stat/run sessions and use GLME to
% compare statistically.

% because different sessions have different directions (both direction and
% #) we treat each direction of motion independently (i.e. a neuron can
% contribute up to 4 tuning curves).


for iarea = 1:numel(areas)
    areaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym], areas(iarea)));
   
    % stat
    allr2 = [areaUnits.r2_stat]';
    allr2p = [areaUnits.r2pval_stat]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),4,1);
    allDir = repmat([1:4]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,4,1);

    statUnitsr2 = cat(1,areaUnits.r2_stat);
    ta(iarea).nStatUnits = sum(any(~isnan(statUnitsr2),2));
    
    % remove instances where insufficient trials to calculate
    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allDir(idx)=[];
    allSigRF(idx)=[];

    ta(iarea).r2_stat = allr2;
    ta(iarea).allr2p_stat = allr2p;
    ta(iarea).areaVec_stat = allArea;
    ta(iarea).sessionVec_stat = allSession;
    ta(iarea).stateVec_stat = zeros(size(ta(iarea).sessionVec_stat));
    ta(iarea).dirVec_stat = allDir;
    ta(iarea).sigRFVec_stat = allSigRF;
    ta(iarea).nTuningCurves_stat = numel(allr2p);
    ta(iarea).nTuned_stat = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_stat = ta(iarea).nTuned_stat/ta(iarea).nTuningCurves_stat;
    ta(iarea).tunedFlag_stat = ta(iarea).r2_stat>r2_thresh & ta(iarea).allr2p_stat<r2p_thresh;


    % run
    allr2 = [areaUnits.r2_run]';
    allr2p = [areaUnits.r2pval_run]';

    allArea = repelem(areas(iarea),numel(allr2),1);
    allSession = repelem(cat(1,areaUnits.sessionID),4,1);
    allDir = repmat([1:4]',numel(areaUnits),1);
    allSigRF = repelem([areaUnits.p_value_rf]'<=0.05,4,1);

    runUnitsr2 = cat(1,areaUnits.r2_run);
    ta(iarea).nRunUnits = sum(any(~isnan(runUnitsr2),2));
   
    % remove instances where insufficient trials to calculate
    idx = find(isnan(allr2));
    allr2(idx) = [];
    allr2p(idx) = [];
    allArea(idx) = [];
    allSession(idx) = [];
    allDir(idx)=[];
    allSigRF(idx)=[];

    ta(iarea).r2_run = allr2;
    ta(iarea).allr2p_run = allr2p;
    ta(iarea).areaVec_run = allArea;
    ta(iarea).sessionVec_run = allSession;
    ta(iarea).stateVec_run = ones(size(ta(iarea).sessionVec_run));
    ta(iarea).dirVec_run = allDir;
    ta(iarea).sigRFVec_run = allSigRF;
    ta(iarea).nTuningCurves_run = numel(allr2p);
    ta(iarea).nTuned_run = sum(allr2>=r2_thresh & allr2p<=r2p_thresh);
    ta(iarea).pTuned_run = ta(iarea).nTuned_run/ta(iarea).nTuningCurves_run;
    ta(iarea).tunedFlag_run = ta(iarea).r2_run>r2_thresh & ta(iarea).allr2p_run<r2p_thresh;

    % combine stat and run into single vectors
    ta(iarea).stateVec = categorical(cat(1,ta(iarea).stateVec_stat,ta(iarea).stateVec_run));
    ta(iarea).tunedVec = cat(1,ta(iarea).tunedFlag_stat,ta(iarea).tunedFlag_run);
    ta(iarea).sessionVec = categorical(cat(1,ta(iarea).sessionVec_stat,ta(iarea).sessionVec_run));
    ta(iarea).dirVec = categorical(cat(1,ta(iarea).dirVec_stat, ta(iarea).dirVec_run));
    ta(iarea).RFVec = categorical(cat(1,ta(iarea).sigRFVec_stat, ta(iarea).sigRFVec_run));
    ta(iarea).areaVec = categorical(cat(1,ta(iarea).areaVec_stat, ta(iarea).areaVec_run));

end

%% plot p(tuned) for each area -  stat, run and run-stat
figure
imagesc([[ta.pTuned_stat]; [ta.pTuned_run]]), colormap(cmocean('tempo')), caxis([0 0.5]), colorbar
ax=gca; ax.XTick=1:8; ax.XTickLabel=areas;
figure
imagesc([[ta.pTuned_run]-[ta.pTuned_stat]]), colormap(crameri('vik')), caxis([-0.25 0.25]), colorbar
ax=gca; ax.XTick=1:8; ax.XTickLabel=areas;

%% create table for glme and perform model selection

stateVec =cat(1,ta.stateVec);
tunedVec = cat(1,ta.tunedVec);
sessionVec = cat(1,ta.sessionVec);
dirVec = cat(1,ta.dirVec);
RFVec = cat(1,ta.RFVec);
areaVec = cat(1,ta.areaVec);

tbl = table(tunedVec, stateVec, sessionVec,dirVec,RFVec,areaVec,...
    'VariableNames', {'tunedFlag', 'state', 'session','dir','RF','area'});

f = 'tunedFlag ~ -1 + area:state + (1|session) + (1|dir) + (1|RF)';
glme_full = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

f = 'tunedFlag ~ -1 + area:state + (1|dir) + (1|RF)';
glme_red = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

f = 'tunedFlag ~ -1 + area:state + (1|RF)';
glme_red2 = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

f = 'tunedFlag ~ -1 + area:state';
glme_red3 = fitglme(tbl,f,'DummyVarCoding','full','Distribution','binomial');

% to compare models: e.g. compare(glme_full,glme_red)

glme = glme_red3

%%  F-tests to test effect of state on p(tuned) in each area

for iarea = 1:8
    H=zeros(1,16);
    H(iarea*2-1)=-1;
    H(iarea*2)=1;
    H;
    ta(iarea).pValue = coefTest(glme,H);
end

% mult-comparisons correction
p_adj = holmbonferroni([ta.pValue]);
format_p = format_p_values(p_adj)


%% plot p(tuned) difference maps between areas with significant tests, separately for stationary and locomotion

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
        statP(iarea1,iarea2) = coefTest(glme,H);

        H=zeros(1,16);
        H(iarea1*2)=-1;
        H(iarea2*2)=1;
        runP(iarea1,iarea2) = coefTest(glme,H);

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

%% tuning strenth distributions

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
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
    nstat(iarea)=numel(statvals);
    nrun(iarea)=numel(runvals);
    title(areas(iarea))
    ylabel('1-CDF'), xlabel('Tuning Strength (R^2)')
end

[nstat; nrun]

%% LME to test effect of state on tuning strength

allAreaUnits = goodUnits(ismember([goodUnits.ecephys_structure_acronym], areas));

stat_r2 = [allAreaUnits.r2_stat];
stat_r2p = [allAreaUnits.r2pval_stat];
stat_area = repelem([allAreaUnits.ecephys_structure_acronym],1,4);
stat_session = repelem([allAreaUnits.sessionID],1,4);
stat_dir = repmat(1:4,1,numel(allAreaUnits));
stat_state = zeros(size(stat_r2));
stat_sigRF = repelem([allAreaUnits.p_value_rf]<=0.05,1,4);

run_r2 = [allAreaUnits.r2_run];
run_r2p = [allAreaUnits.r2pval_run];
run_area = repelem([allAreaUnits.ecephys_structure_acronym],1,4);
run_session = repelem([allAreaUnits.sessionID],1,4);
run_dir = repmat(1:4,1,numel(allAreaUnits));
run_state = ones(size(run_r2));
run_sigRF = repelem([allAreaUnits.p_value_rf]<=0.05,1,4);

all_r2 = [stat_r2, run_r2]; all_r2(all_r2<0)=0; % set to zero floor
all_r2p = [stat_r2p, run_r2p];
all_area = [stat_area, run_area];
all_session = [stat_session, run_session];
all_dir = [stat_dir, run_dir];
all_state = [stat_state, run_state];
tunedFlag = all_r2>=r2_thresh & all_r2p<=r2p_thresh;
all_sigRF = [stat_sigRF, run_sigRF];

idx = find(isnan(all_r2));
all_r2(idx) = [];
all_r2p(idx) = [];
all_area(idx) = [];
all_session(idx) = [];
all_dir(idx) = [];
all_state(idx) = [];
tunedFlag(idx) = [];
all_sigRF(idx) = [];



tbl = table(tunedFlag', all_r2(:), all_area', categorical(all_session'),...
    categorical(all_dir'), categorical(all_state'), categorical(all_sigRF(:)),...
    'VariableNames', {'tunedFlag', 'tuningStrength', 'area', 'session', 'dir', 'state', 'sigRF'});

f = 'tuningStrength ~ -1 + area:state + (1|session) + (1|dir) + (1|sigRF)';
lme_full = fitlme(tbl,f,'DummyVarCoding','full');

% f = 'tuningStrength ~ -1 + area:state +  (1|dir) + (1|sigRF)';
% lme_red = fitlme(tbl,f,'DummyVarCoding','full');
% 
% f = 'tuningStrength ~ -1 + area:state +  (1|sigRF)';
% lme_red2 = fitlme(tbl,f,'DummyVarCoding','full');
% 
% f = 'tuningStrength ~ -1 + area:state';
% lme_red3 = fitlme(tbl,f,'DummyVarCoding','full');

lme = lme_full;

areaOrder = [3, 7, 5, 8, 1, 2, 4, 6];
[~, sortidx] = sort(areaOrder);
H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1) = -1; H_temp(iarea1+8) = 1;
    pvals(iarea1) = coefTest(lme,H_temp);
end

pvals = pvals(sortidx);
pvals_adj = holmbonferroni(pvals)

%% lme for tuning strength, only for tuned cells 

tbl = table(tunedFlag', all_r2(:), all_area', categorical(all_session'),...
    categorical(all_dir'), categorical(all_state'), categorical(all_sigRF(:)),...
    'VariableNames', {'tunedFlag', 'tuningStrength', 'area', 'session', 'dir', 'state', 'sigRF'});

tbl_tuned = tbl(tbl.tunedFlag,:);
f = 'tuningStrength ~ -1 + area:state + (1|session) + (1|dir) + (1|sigRF)';
lme = fitlme(tbl_tuned,f,'DummyVarCoding','full');

areaOrder = [3, 7, 5, 8, 1, 2, 4, 6];
[~, sortidx] = sort(areaOrder);
H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1) = -1; H_temp(iarea1+8) = 1;
    pvals(iarea1) = coefTest(lme,H_temp);
end

pvals = pvals(sortidx);
pvals_adj = holmbonferroni(pvals)

% get actual tuning strength values
for iarea = 1:8
    temp_tbl = tbl_tuned(strcmp(tbl_tuned.area, areas(iarea)),:);
    statta(iarea).tuned_r2 = temp_tbl.tuningStrength(double(string(temp_tbl.state))==0);
    runta(iarea).tuned_r2 = temp_tbl.tuningStrength(double(string(temp_tbl.state))==1);
end

% mean+-SEM for stat and run 
tuningStrengthForTunedCells = [cellfun(@mean, {statta.tuned_r2}); cellfun(@sem,  {statta.tuned_r2});...
cellfun(@mean, {runta.tuned_r2}); cellfun(@sem,  {runta.tuned_r2})]



%% Supp. Figure 5 plot tuning strength cdfs, all areas, by state

figure
    subplot(121), 
    hold on
    for iarea = 1:8
    statvals = ta(iarea).r2_stat; statvals(statvals<0)=-.1;
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;
    [f, x] = ecdf(statvals);
    plot(x,1-f,'Color', areacols(iarea,:).*0.7);
    xlim([0 1])
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
    end
    title('stat')
    ylabel('1-CDF')
    xlabel('Tuning Strength (R^2)')

        subplot(122), 
    hold on
    for iarea = 1:8
    runvals = ta(iarea).r2_run; runvals(runvals<0)=-.1;
    [f, x] = ecdf(runvals);
    plot(x,1-f,'Color', areacols(iarea,:));
    xlim([0 1])
    ylim([0 0.6])
    
    ax=gca; ax.YTick = 0:0.1:0.6; ax.XTick = 0:0.1:1;
    grid on
    defaultAxesProperties(gca,false)
            title('run')

    end

