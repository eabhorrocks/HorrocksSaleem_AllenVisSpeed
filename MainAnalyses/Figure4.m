%% Figure 4 - Decoding analysis

%% load the data

dataDir = 'D:\AllenMatFiles\FC';

dataFiles = dir(fullfile(dataDir, 'dmPID*'));
for ifile = 1:numel(dataFiles)
    dataFiles(ifile).sessionID = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'));
    dataFiles(ifile).filename = fullfile(dataFiles(ifile).folder, dataFiles(ifile).name);
end

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
allColors = tab20(20);
areacols = allColors([1 3 5 7 9 11 13 17],:);



for isession = 1:numel(dataFiles)

    load(dataFiles(isession).filename)
    session(isession).stat = stat;
    session(isession).run=run;
end

%% format data for plotting and stats
% each pCorrect with area, session and dir

nDirs=4;

for iarea=1:8
ta(iarea).statMeanPerf=[];
ta(iarea).statDir = [];
ta(iarea).statSesh = [];

ta(iarea).runMeanPerf=[];
ta(iarea).runDir = [];
ta(iarea).runSesh = [];
end

% stat [erf
for isession = 1:numel(dataFiles)
    for iarea = 1:8
            for iperm = 1:numel(session(isession).stat.ta(iarea).perm)
                ta(iarea).statMeanPerf = cat(2,ta(iarea).statMeanPerf, session(isession).stat.ta(iarea).perm(iperm).dir.mean_pCorrect);
                ta(iarea).statDir = cat(2,ta(iarea).statDir, 1:nDirs);
                ta(iarea).statSesh = cat(2,ta(iarea).statSesh, repelem(isession,1,nDirs));
            end
    end
end

% run perf
for isession = 1:numel(dataFiles)
    for iarea = 1:8
            for iperm = 1:numel(session(isession).run.ta(iarea).perm)
                ta(iarea).runMeanPerf = cat(2,ta(iarea).runMeanPerf, session(isession).run.ta(iarea).perm(iperm).dir.mean_pCorrect);
                ta(iarea).runDir = cat(2,ta(iarea).runDir, 1:nDirs);
                ta(iarea).runSesh = cat(2,ta(iarea).runSesh, repelem(isession,1,nDirs));
            end
    end
end

%% paired box plot
figure, hold on
nestedCellArray = {{ta.statMeanPerf},{ta.runMeanPerf}};
colorCellArray = {{areacols.*0.7}, {areacols}};
subGroupSpacing = [-0.15, +0.15];
[bp, a] = pairedBoxPlot(nestedCellArray, subGroupSpacing, colorCellArray);
ax= gca;
ax.XTick = 1:8;
ylim([0 0.7])
yline(1/7,'k:') % chance

%% LME statistical analysis

stat_perf = [ta.statMeanPerf];
stat_dir = [ta.statDir];
stat_session = [ta.statSesh];
stat_area = repelem(areas, cellfun(@numel, {ta.statMeanPerf}));
stat_state = zeros(size(stat_perf));

run_perf = [ta.runMeanPerf];
run_dir = [ta.runDir];
run_session = [ta.runSesh];
run_area = repelem(areas, cellfun(@numel, {ta.runMeanPerf}));
run_state = ones(size(run_perf));

perfVec = [stat_perf(:); run_perf(:)];
stateVec = [stat_state(:); run_state(:)];
sessionVec = [stat_session(:); run_session(:)];
areaVec = [stat_area(:); run_area(:)];
dirVec = [stat_dir(:); run_dir(:)];

tbl = table(perfVec, stateVec, sessionVec, areaVec, dirVec,...
    'VariableNames', {'perf', 'state', 'session', 'area', 'dir'});

tbl.session = categorical(tbl.session);
tbl.state = categorical(tbl.state);
tbl.dir = categorical(tbl.dir);

f = 'perf ~ -1 + area:state + (1|dir)';
lme = fitlme(tbl, f, 'DummyVarCoding', 'full')

% get effect of state for each area
H = zeros(1,16);
pvals = nan(1,8);
for iarea1 = 1:8
    H_temp = H;
    H_temp(iarea1*2-1) = -1; H_temp(iarea1*2) = 1;
    pvals(iarea1) = coefTest(lme,H_temp); %
end


pvals_adj = holmbonferroni(pvals)

%% report mean+-SEM values for decoding performance

[cellfun(@(x) nanmean(x,2), {ta.statMeanPerf});...
    cellfun(@(x) nansem(x,2),  {ta.statMeanPerf});...
    cellfun(@(x) nanmean(x,2), {ta.runMeanPerf});...
    cellfun(@(x) nansem(x,2),  {ta.runMeanPerf})]



