function allen_analyseSesssion_DotMotionSpeedTuning(inputFileName, outputFileName)

%% load session
load(inputFileName)

%% options
statMeanThresh = 0.5; %0.5
statUpperThresh = 3;
runMeanThresh = 3;
runLowerThresh = 0.5;
reqTrials = 8;

doPSTH=true;

% binned spike count options
scopts.intervalStart = 0;
scopts.intervalEnd = 1;
scopts.binSpacing=1;
scopts.uniqueVals = {[-45 0 45 90], [0, 0.0005, 0.0010, 0.0050, 0.0100, 0.0200, 0.0400]};
alldirs = [-45 0 45 90];
speeds = [0, 0.0005, 0.0010, 0.0050, 0.0100, 0.0200, 0.0400];

blscopts = scopts;
blscopts.intervalStart = -0.2;
blscopts.intervalEnd = 0;
blscopts.binSpacing=0.2;
blscopts.uniqueVals = {[-45 0 45 90], [0, 0.0005, 0.0010, 0.0050, 0.0100, 0.0200, 0.0400]};

% r2 calc options
kfold = 3; nPerms = 10; randFlag = 1; validMeans = true; nShuffle = 100;

% get valid trials
dmtrials = trials(strcmp([trials.stimulus_name], 'dot_motion'));
dmtrials_valid = dmtrials(~[dmtrials.containsInvalidTime]);

% V1 dynamics
run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials_valid.runTrace}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials_valid.runTrace}));

% simpler
%run_idx = find(cellfun(@(x) prop(x>0.5)>=0 & mean(x)>3, {dmtrials_valid.runTrace}));
%stat_idx = find(cellfun(@(x) prop(x<3)>=0 & mean(x)<0.5, {dmtrials_valid.runTrace}));

statTrials = dmtrials_valid(stat_idx);
runTrials = dmtrials_valid(run_idx);


% MI options
optimiseBins = false; nMCSamples = 1000; sigFlag = false; correction = 'MLE'; nBinLim = 7;

[units.spikeCounts_stat] = deal(cell(4,7));
[units.spikeCounts_run] = deal(cell(4,7));
[units.baselineSpikeCounts_stat] = deal(cell(4,7));
[units.baselineSpikeCounts_run] = deal(cell(4,7));
[units.tuning_stat] = deal(nan(4,7));
[units.tuning_run] = deal(nan(4,7));
[units.dynamicRange_stat] = deal(nan(1,4));
[units.dynamicRange_run] = deal(nan(1,4));
[units.fanoFactor_stat] = deal(nan(1,4));
[units.fanoFactor_run] = deal(nan(1,4));


%% get binned spike counts

if numel(unique([dmtrials.Speed]==7)) % 2 sessions don't have the 7 speeds

    if ~isempty(statTrials)
        % get binned spike counts, tuning, dynamic range and fano factor
        [anUnits, cond] = getBinnedSpikeCounts(statTrials, units, {'Dir', 'Speed'}, scopts);
        for iunit = 1:numel(units)
            units(iunit).spikeCounts_stat = anUnits(iunit).allSpikes;
            units(iunit).spikeCounts_stat = cellfun(@(x) x', units(iunit).spikeCounts_stat,'UniformOutput',false);
            units(iunit).tuning_stat = cellfun(@mean, units(iunit).spikeCounts_stat);
            units(iunit).dynamicRange_stat = range(units(iunit).tuning_stat,2);
            units(iunit).fanoFactor_stat = mean(cellfun(@(x) var(x)/mean(x), units(iunit).spikeCounts_stat),2);
        end
    end

    if ~isempty(runTrials)
        [anUnits, cond] = getBinnedSpikeCounts(runTrials, units, {'Dir', 'Speed'}, scopts);
        for iunit = 1:numel(units)
            units(iunit).spikeCounts_run = anUnits(iunit).allSpikes;
            units(iunit).spikeCounts_run = cellfun(@(x) x', units(iunit).spikeCounts_run,'UniformOutput',false);
            units(iunit).tuning_run = cellfun(@mean, units(iunit).spikeCounts_run);
            units(iunit).dynamicRange_run = range(units(iunit).tuning_run,2);
            units(iunit).fanoFactor_run = mean(cellfun(@(x) var(x)/mean(x), units(iunit).spikeCounts_run),2);
        end
    end


    %% get baseline firing rate estimates
    % binned spike count options


    if ~isempty(statTrials)
        % get binned spike counts, tuning, dynamic range and fano factor
        [anUnits, cond] = getBinnedSpikeCounts(statTrials, units, {'Dir', 'Speed'}, blscopts);
        for iunit = 1:numel(units)
            units(iunit).baselineSpikeCounts_stat = anUnits(iunit).allSpikes;
            units(iunit).baselineSpikeCounts_stat = cellfun(@(x) x', units(iunit).baselineSpikeCounts_stat,'UniformOutput',false);
        end
    end

    if ~isempty(runTrials)
        [anUnits, cond] = getBinnedSpikeCounts(runTrials, units, {'Dir', 'Speed'}, blscopts);
        for iunit = 1:numel(units)
            units(iunit).baselineSpikeCounts_run = anUnits(iunit).allSpikes;
            units(iunit).baselineSpikeCounts_run = cellfun(@(x) x', units(iunit).baselineSpikeCounts_run,'UniformOutput',false);
        end
    end

    %% set nans as default values



    [units.r2_stat] = deal(nan(1,4));
    [units.r2_run] = deal(nan(1,4));
    [units.r2pval_stat] = deal(nan(1,4));
    [units.r2pval_run] = deal(nan(1,4));

    [units.MI_stat] = deal(nan(1,4));
    [units.SSI_stat] = deal(nan(4,7));
    [units.MI_run] = deal(nan(1,4));
    [units.SSI_run] = deal(nan(4,7));

    [units.gaussParams_stat] = deal(nan(4));
    [units.gaussParams_run] = deal(nan(4));
    [units.gaussChar_stat] = deal(nan(1,4));
    [units.gaussChar_run] = deal(nan(1,4));
    [units.gaussR2_stat] = deal(nan(1,4));
    [units.gaussR2_run] = deal(nan(1,4));
    [units.prefSpeed_stat] = deal(nan(1,4));
    [units.prefSpeed_run] = deal(nan(1,4));

    for iunit = 1:numel(units)
        for idir = 1:4
            units(iunit).statdirpsth(idir).PSTH = [];
            units(iunit).statdirpsth(idir).info = [];
            units(iunit).rundirpsth(idir).PSTH = [];
            units(iunit).rundirpsth(idir).info = [];
        end
    end



    %% do the analysis


    %% stationary trials
    if ~isempty(statTrials)

        for idir = 1:4
            % get trials
            % if statTrials
            if min(cellfun(@numel, units(1).spikeCounts_stat(idir,:)))>=reqTrials

                %
                for iunit=1:numel(units)

                    % downsample trials and get tuning strength
                    % spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    %     units(iunit).spikeCounts_stat(idir,:),'UniformOutput',false);

                    spca = cellfun(@(x) x(1:reqTrials),...
                        units(iunit).spikeCounts_stat(idir,:),'UniformOutput',false);

                    [units(iunit).r2_stat(idir), units(iunit).r2pval_stat(idir)] = ...
                        calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

                    % info theory with downsampled trials
                    [units(iunit).MI_stat(1,idir), ~, ~, units(iunit).SSI_stat(idir,:)] =...
                        calcMI(spca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);

                    % gaussian fit and preferred speed with all trials
                    [units(iunit).gaussParams_stat(idir,:), units(iunit).gaussChar_stat(idir), units(iunit).gaussR2_stat(idir)] = ...
                        fitGaussianTemplates_tuning(units(iunit).spikeCounts_stat(idir,:),0.5,false);

                    if units(iunit).gaussChar_stat(1,idir) == 4 % if inverted, pref speed is min fr.
                        [~, units(iunit).prefSpeed_stat(1,idir)] = min(units(iunit).tuning_stat(idir,:));
                    else % max fr
                        [~, units(iunit).prefSpeed_stat(1,idir)] = max(units(iunit).tuning_stat(idir,:));
                    end
                end





                %% PSTHs

                if doPSTH

                    areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

                    % options for PSTH function
                    binWidth = 10; % ms
                    options.smoothType = 'gaussian';
                    options.smoothWidth = 175; % in ms
                    options.preTime = 500; % time in ms before stim onset
                    options.postTime = 500; % time in ms after stim onset
                    options.plot = false; % don't make plots if processing all units
                    % psth reliability options
                    options.getReliability = false;
                    options.minReliTrials = 9; % min number of trials to do reliability
                    options.nReliPerms = 50; % number of iterations to compute reliability
                    options.nShuffle = 20; % number of shuffles per cv iteration
                    options.maxStretch = 20; % max sample value for dtw function (in ms, must be divisible by binWidth).

                    % stationary trials

                    clear sints
                    for ispeed = 1:numel(speeds)
                        t = statTrials([statTrials.Speed]==speeds(ispeed) & [statTrials.Dir]==alldirs(idir)); % get the relevant trials for this condition
                        sints{ispeed}= [vertcat(t.start_time), vertcat(t.stop_time)].*1000; % stimulus in interval times, convert to ms
                    end


                    for iunit = 1:numel(units)
                        if ismember(units(iunit).ecephys_structure_acronym, areas)

                            [units(iunit).statdirpsth(idir).PSTH, units(iunit).statdirpsth(idir).info] =...
                                spikeTimesPSTH_paper(units(iunit).spiketimes.*1000,sints,binWidth,options); % convert spike times to ms
                        end
                    end


                    % fit gaussian descriptive functions to characterise onset and offset features

                    earlyidx = 51:80; % time idx of onset period
                    offsetidx = 151:200; % offset period
                    customOffset_onset = 5; % a parameter related to the mean of the gaussian
                    customOffset_offset = 5;

                    psth_len = 200;  % element length of PSTH

                    for iunit = 1:numel(units)
                        if ismember(units(iunit).ecephys_structure_acronym, areas)



                            for ispeed = 1:7
                                psth_norm = normalize(units(iunit).statdirpsth(idir).PSTH(ispeed).psth,'range'); % normalise PSTH
                                while length(psth_norm)<psth_len
                                    psth_norm(end+1) = psth_norm(end);
                                end
                                onsetPSTH = num2cell(psth_norm(earlyidx));
                                offsetPSTH = num2cell(psth_norm(offsetidx));

                                [~, units(iunit).statdirpsth(idir).info(ispeed).onset_char, ...
                                    units(iunit).statdirpsth(idir).info(ispeed).onset_charR2] =...
                                    fitGaussianTemplates_psth(onsetPSTH,customOffset_onset,false);

                                [~, units(iunit).statdirpsth(idir).info(ispeed).offset_char,...
                                    units(iunit).statdirpsth(idir).info(ispeed).offset_charR2] =...
                                    fitGaussianTemplates_psth(offsetPSTH,customOffset_offset,false);

                            end
                        end

                    end
                end

            end

        end
    end



    %% Running trials
    if ~isempty(runTrials)

        for idir = 1:4
            % get trials
            % if statTrials
            if min(cellfun(@numel, units(1).spikeCounts_run(idir,:)))>=reqTrials

                %
                for iunit=1:numel(units)

                    % downsample trials and get tuning strength
                    % spca = cellfun(@(x) x(randsample(1:numel(x),reqTrials)),...
                    %     units(iunit).spikeCounts_run(idir,:),'UniformOutput',false);

                    spca = cellfun(@(x) x(1:reqTrials),...
                        units(iunit).spikeCounts_run(idir,:),'UniformOutput',false);

                    [units(iunit).r2_run(idir), units(iunit).r2pval_run(idir)] = ...
                        calc_kfold_R2(spca, kfold, nPerms, randFlag, validMeans, nShuffle);

                    % info theory with downsampled trials
                    [units(iunit).MI_run(1,idir), ~, ~, units(iunit).SSI_run(idir,:)] =...
                        calcMI(spca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);

                    % gaussian fit and preferred speed with all trials
                    [units(iunit).gaussParams_run(idir,:), units(iunit).gaussChar_run(idir), units(iunit).gaussR2_run(idir)] = ...
                        fitGaussianTemplates_tuning(units(iunit).spikeCounts_run(idir,:),0.5,false);

                    if units(iunit).gaussChar_run(1,idir) == 4 % if inverted, pref speed is min fr.
                        [~, units(iunit).prefSpeed_run(1,idir)] = min(units(iunit).tuning_run(idir,:));
                    else % max fr
                        [~, units(iunit).prefSpeed_run(1,idir)] = max(units(iunit).tuning_run(idir,:));
                    end
                end


                %% PSTHs

                if doPSTH

                    areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};



                    % options for PSTH function
                    binWidth = 10; % ms
                    options.smoothType = 'gaussian';
                    options.smoothWidth = 175; % in ms
                    options.preTime = 500; % time in ms before stim onset
                    options.postTime = 500; % time in ms after stim onset
                    options.plot = false; % don't make plots if processing all units
                    % psth reliability options
                    options.getReliability = true;
                    options.minReliTrials = 9; % min number of trials to do reliability
                    options.nReliPerms = 50; % number of iterations to compute reliability
                    options.nShuffle = 20; % number of shuffles per cv iteration
                    options.maxStretch = 20; % max sample value for dtw function (in ms, must be divisible by binWidth).

                    % runionary trials

                    clear sints
                    for ispeed = 1:numel(speeds)
                        t = runTrials([runTrials.Speed]==speeds(ispeed) & [runTrials.Dir]==alldirs(idir)); % get the relevant trials for this condition
                        sints{ispeed}= [vertcat(t.start_time), vertcat(t.stop_time)].*1000; % stimulus in interval times, convert to ms
                    end


                    for iunit = 1:numel(units)
                        if ismember(units(iunit).ecephys_structure_acronym, areas)

                            [units(iunit).rundirpsth(idir).PSTH, units(iunit).rundirpsth(idir).info] =...
                                spikeTimesPSTH_paper(units(iunit).spiketimes.*1000,sints,binWidth,options); % convert spike times to ms
                        end
                    end



                    % fit gaussian descriptive functions to characterise onset and offset features

                    earlyidx = 51:80; % time idx of onset period
                    offsetidx = 151:200; % offset period
                    customOffset_onset = 5; % a parameter related to the mean of the gaussian
                    customOffset_offset = 5;

                    psth_len = 200;  % element length of PSTH

                    for iunit = 1:numel(units)
                        if ismember(units(iunit).ecephys_structure_acronym, areas)

                            for ispeed = 1:7
                                psth_norm = normalize(units(iunit).rundirpsth(idir).PSTH(ispeed).psth,'range'); % normalise PSTH
                                while length(psth_norm)<psth_len
                                    psth_norm(end+1) = psth_norm(end);
                                end
                                onsetPSTH = num2cell(psth_norm(earlyidx));
                                offsetPSTH = num2cell(psth_norm(offsetidx));

                                [~, units(iunit).rundirpsth(idir).info(ispeed).onset_char, ...
                                    units(iunit).rundirpsth(idir).info(ispeed).onset_charR2] =...
                                    fitGaussianTemplates_psth(onsetPSTH,customOffset_onset,false);

                                [~, units(iunit).rundirpsth(idir).info(ispeed).offset_char,...
                                    units(iunit).rundirpsth(idir).info(ispeed).offset_charR2] =...
                                    fitGaussianTemplates_psth(offsetPSTH,customOffset_offset,false);

                            end
                        end

                    end
                end

            end

        end
    end

    %% save data

    save(outputFileName, 'units', 'trials', 'runInfo', 'pupilInfo', 'sessionInfo')


end


end