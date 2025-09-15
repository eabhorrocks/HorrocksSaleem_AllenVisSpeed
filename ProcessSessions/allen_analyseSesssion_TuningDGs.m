function allen_analyseSesssion_TuningDGs(inputFileName, outputFileName)

%% load session
load(inputFileName)

%% options
statMeanThresh = 1; %0.5
statUpperThresh = 3;
runMeanThresh = 3;
runLowerThresh = 0.5;
reqTrials = 8;

doPSTH=false;


%% get valid trials
dmtrials = trials(strcmp([trials.stimulus_name], 'drifting_gratings'));
dmtrials_valid = dmtrials(~[dmtrials.containsInvalidTime]);
dmtrials_valid = dmtrials(~isnan([dmtrials.orientation]));



% if numel(unique([dmtrials.Speed]==7)) % 2 sessions don't have the 7 speeds

run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials_valid.runTrace}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials_valid.runTrace}));

statTrials = dmtrials_valid(stat_idx);
runTrials = dmtrials_valid(run_idx);


uniqueDirs = unique([dmtrials_valid.orientation]);
nDirs=numel(uniqueDirs);
uniqueTF = unique([dmtrials_valid.temporal_frequency]);
nTF = numel(uniqueTF);



[units.spikeCounts_stat] = deal(cell(4,7));
[units.spikeCounts_run] = deal(cell(4,7));
[units.tuning_stat] = deal(nan(4,7));
[units.tuning_run] = deal(nan(4,7));
[units.dynamicRange_stat] = deal(nan(1,4));
[units.dynamicRange_run] = deal(nan(1,4));
[units.fanoFactor_stat] = deal(nan(1,4));
[units.fanoFactor_run] = deal(nan(1,4));

%% get binned spike counts

% binned spike count options
scopts.intervalStart = 0;
scopts.intervalEnd = 2; % 2s trials
scopts.binSpacing=2;
scopts.uniqueVals = {uniqueDirs, uniqueTF};

% r2 calc options
kfold = 3; nPerms = 10; randFlag = 1; validMeans = true; nShuffle = 100;



if ~isempty(statTrials)
    % get binned spike counts, tuning, dynamic range and fano factor
    [anUnits, cond] = getBinnedSpikeCounts(statTrials, units, {'orientation', 'temporal_frequency'}, scopts);
    for iunit = 1:numel(units)
        units(iunit).spikeCounts_stat = anUnits(iunit).allSpikes;
        units(iunit).spikeCounts_stat = cellfun(@(x) x', units(iunit).spikeCounts_stat,'UniformOutput',false);
        units(iunit).tuning_stat = cellfun(@mean, units(iunit).spikeCounts_stat);
        units(iunit).dynamicRange_stat = range(units(iunit).tuning_stat,2);
        units(iunit).fanoFactor_stat = mean(cellfun(@(x) var(x)/mean(x), units(iunit).spikeCounts_stat),2);
    end
end

if ~isempty(runTrials)
    [anUnits, cond] = getBinnedSpikeCounts(runTrials, units, {'orientation', 'temporal_frequency'}, scopts);
    for iunit = 1:numel(units)
        units(iunit).spikeCounts_run = anUnits(iunit).allSpikes;
        units(iunit).spikeCounts_run = cellfun(@(x) x', units(iunit).spikeCounts_run,'UniformOutput',false);
        units(iunit).tuning_run = cellfun(@mean, units(iunit).spikeCounts_run);
        units(iunit).dynamicRange_run = range(units(iunit).tuning_run,2);
        units(iunit).fanoFactor_run = mean(cellfun(@(x) var(x)/mean(x), units(iunit).spikeCounts_run),2);
    end
end


%% set nans as default values



[units.r2dir_stat] = deal(nan(1,nTF));
[units.r2dir_run] = deal(nan(1,nTF));
[units.r2dirpval_stat] = deal(nan(1,nTF));
[units.r2dirpval_run] = deal(nan(1,nTF));

[units.r2tf_stat] = deal(nan(1,nDirs));
[units.r2tf_run] = deal(nan(1,nDirs));
[units.r2tfpval_stat] = deal(nan(1,nDirs));
[units.r2tfpval_run] = deal(nan(1,nDirs));

% [units.MI_stat] = deal(nan(1,4));
% [units.SSI_stat] = deal(nan(4,7));
% [units.MI_run] = deal(nan(1,4));
% [units.SSI_run] = deal(nan(4,7));
%
% [units.gaussParams_stat] = deal(nan(4));
% [units.gaussParams_run] = deal(nan(4));
% [units.gaussChar_stat] = deal(nan(1,4));
% [units.gaussChar_run] = deal(nan(1,4));
% [units.gaussR2_stat] = deal(nan(1,4));
% [units.gaussR2_run] = deal(nan(1,4));
% [units.prefSpeed_stat] = deal(nan(1,4));
% [units.prefSpeed_run] = deal(nan(1,4));




%% do the analysis

%% TF tuning
if ~isempty(statTrials)

    for idir = 1:nDirs % tf tuning for each dir

        if min(cellfun(@numel, units(1).spikeCounts_stat(idir,:)))>=reqTrials

            %
            for iunit=1:numel(units)

                % downsample trials and get tuning strength
                spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    units(iunit).spikeCounts_stat(idir,:),'UniformOutput',false);

                [units(iunit).r2tf_stat(idir), units(iunit).r2tfpval_stat(idir)] = ...
                    calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

            end
        end
    end
end

if ~isempty(runTrials)

    for idir = 1:nDirs % tf tuning for each dir

        if min(cellfun(@numel, units(1).spikeCounts_run(idir,:)))>=reqTrials

            %
            for iunit=1:numel(units)

                % downsample trials and get tuning strength
                spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    units(iunit).spikeCounts_run(idir,:),'UniformOutput',false);

                [units(iunit).r2tf_run(idir), units(iunit).r2tfpval_run(idir)] = ...
                    calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

            end
        end
    end
end


%% Dir tuning

if ~isempty(statTrials)

    for itf = 1:nTF % dir tuning for each tf
        if min(cellfun(@numel, units(1).spikeCounts_stat(:,itf)))>=reqTrials

            %
            for iunit=1:numel(units)

                % downsample trials and get tuning strength
                spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    units(iunit).spikeCounts_stat(:,itf),'UniformOutput',false);

                [units(iunit).r2dir_stat(itf), units(iunit).r2dirpval_stat(itf)] = ...
                    calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

            end
        end
    end
end


if ~isempty(runTrials)

    for itf = 1:nTF % dir tuning for each tf
        if min(cellfun(@numel, units(1).spikeCounts_run(:,itf)))>=reqTrials

            %
            for iunit=1:numel(units)

                % downsample trials and get tuning strength
                spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    units(iunit).spikeCounts_run(:,itf),'UniformOutput',false);

                [units(iunit).r2dir_run(itf), units(iunit).r2dirpval_run(itf)] = ...
                    calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

            end
        end
    end
end



%% save data

save(outputFileName, 'units', 'trials', 'runInfo', 'pupilInfo', 'sessionInfo')


end


