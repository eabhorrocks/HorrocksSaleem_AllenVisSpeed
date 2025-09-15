function allen_analyseSesssion_PIDSpeedTuning(inputFileName, outputFileName)

%% load session
load(inputFileName)


%% options and vars

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

reqTrials = 8;
reqUnits=10;

% binned spike count options
scopts.intervalStart = 0;
scopts.intervalEnd = 1;
scopts.binSpacing=1;
scopts.uniqueVals = {[-45 0 45 90], [0, 0.0005, 0.0010, 0.0050, 0.0100, 0.0200, 0.0400]};
alldirs = [-45 0 45 90];
speeds = [0, 0.0005, 0.0010, 0.0050, 0.0100, 0.0200, 0.0400];
nConds = numel(speeds);

% get valid trials
dmtrials = trials(strcmp([trials.stimulus_name], 'dot_motion'));
dmtrials_valid = dmtrials(~[dmtrials.containsInvalidTime]);


run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials_valid.runTrace}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials_valid.runTrace}));

statTrials = dmtrials_valid(stat_idx);
runTrials = dmtrials_valid(run_idx);

%% get binned spike counts
if numel(unique([dmtrials.Speed]==7)) % 2 sessions don't have the 7 speeds

    if ~isempty(statTrials)
        % get binned spike counts, tuning, dynamic range and fano factor
        [anUnits, cond] = getBinnedSpikeCounts(statTrials, units, {'Dir', 'Speed'}, scopts);
        for iunit = 1:numel(units)
            units(iunit).spikeCounts_stat = anUnits(iunit).allSpikes;
            units(iunit).spikeCounts_stat = cellfun(@(x) x', units(iunit).spikeCounts_stat,'UniformOutput',false);
        end
    end

    if ~isempty(runTrials)
        [anUnits, cond] = getBinnedSpikeCounts(runTrials, units, {'Dir', 'Speed'}, scopts);
        for iunit = 1:numel(units)
            units(iunit).spikeCounts_run = anUnits(iunit).allSpikes;
            units(iunit).spikeCounts_run = cellfun(@(x) x', units(iunit).spikeCounts_run,'UniformOutput',false);
        end
    end

    clear cond


    %% use good units
    goodUnits = units([units.isi_violations]<=0.1...
        & [units.amplitude_cutoff]<=0.1 & [units.waveform_amplitude]>=50);




    %% Stationary

    scName = 'spikeCounts_stat';
    if ~isempty(statTrials)

        for idir = 1:4

            % check trials available
            minTrials = min(cellfun(@numel, units(1).(scName)(idir,:)));
            if minTrials<reqTrials % if not enough trials, set to nan
                for iarea = 1:8
                    tempAreaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym],areas(iarea)));
                    nUnits = numel(tempAreaUnits);
                    nSubPops = floor(nUnits/reqUnits); % number of subpops we can gert
                    if nSubPops<1, nSubPops=1; end
                    for isubpop = 1:nSubPops
                    ta(iarea).perm(isubpop).dir(idir).mean_pCorrect = nan(1,1);
                    ta(iarea).perm(isubpop).dir(idir).mean_Error = nan(1,1);
                    ta(iarea).perm(isubpop).dir(idir).mean_bias = nan(1,1);
                    end
                end
            else % if enough trials, do decoding

                for iarea = 1:numel(areas)

                    tempAreaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym],areas(iarea)));
                    nUnits = numel(tempAreaUnits);
                    nSubPops = floor(nUnits/reqUnits); % number of subpops we can gert

                    if nSubPops<1 % if not enough units to do any decoding, set to nan
                        for idir2=1:4
                        ta(iarea).perm(1).dir(idir2).mean_pCorrect = nan(1,1);
                        ta(iarea).perm(1).dir(idir2).mean_Error = nan(1,1);
                        ta(iarea).perm(1).dir(idir2).mean_bias = nan(1,1);
                        end
                    else % if enough units, do decoding

                        tempAreaUnits = tempAreaUnits(randperm(nUnits)); % shuffle order of units for the sake of it

                        for isubpop = 1:nSubPops
                            tempUnits = tempAreaUnits(isubpop*reqUnits-(reqUnits-1):(isubpop*reqUnits)); % units for this subpop

                            % get the trials to be used for decoding (downsampling to nTrials)
                            for icond = 1:nConds
                                nUseableTrials = size(tempUnits(1).(scName){idir,icond},1);
                                    randSampTrials = randsample(nUseableTrials,reqTrials);
                                cond{icond}.trials = randSampTrials;
                            end


                            % generate input struct for decoding function
                            units2Decode = []; units2Decode(reqUnits).spikeCell = [];

                            for iunit = 1:numel(tempUnits)
                                for icond = 1:nConds
                                    units2Decode(iunit).spikeCell{icond} =...
                                        tempUnits(iunit).(scName){idir,icond}...
                                        (cond{icond}.trials);
                                end
                            end

                            % do the LOO-CV decoding
                            tuningAxes = {[1:7]};
                            ta(iarea).perm(isubpop).dir(idir).cond = PID_DecodeFromSpikeCellArray(units2Decode, tuningAxes, 0);

                            % get pCorrect and bias of predictions for each speed, for
                            % this subpop
                            for icond = 1:nConds
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).pCorrect = prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds==icond);
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).meanError = mean(abs(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds-icond));
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).bias = prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds>icond) - prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds<icond);
                            end

                            ta(iarea).perm(isubpop).dir(idir).mean_pCorrect = mean([ta(iarea).perm(isubpop).dir(idir).cond.pCorrect]);
                            ta(iarea).perm(isubpop).dir(idir).mean_Error = mean([ta(iarea).perm(isubpop).dir(idir).cond.meanError]);
                            ta(iarea).perm(isubpop).dir(idir).mean_bias = mean([ta(iarea).perm(isubpop).dir(idir).cond.bias]);
                        end

                    end

                end

            end

        end

    else % set all nans if no trials
        for iarea = 1:8
            for idir = 1:4
                ta(iarea).perm(1).dir(idir).mean_pCorrect = nan(1,1);
                ta(iarea).perm(1).dir(idir).mean_Error = nan(1,1);
                ta(iarea).perm(1).dir(idir).mean_bias = nan(1,1);
            end
        end



    end


    stat.ta=ta;
    clear ta;




    %% Locomotion

    scName = 'spikeCounts_run';
    if ~isempty(runTrials)

        for idir = 1:4

            % check trials available
            minTrials = min(cellfun(@numel, units(1).(scName)(idir,:)));
            if minTrials<reqTrials % if not enough trials, set to nan
                for iarea = 1:8
                    tempAreaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym],areas(iarea)));
                    nUnits = numel(tempAreaUnits);
                    nSubPops = floor(nUnits/reqUnits); % number of subpops we can gert
                    if nSubPops<1, nSubPops=1; end
                    for isubpop = 1:nSubPops
                    ta(iarea).perm(isubpop).dir(idir).mean_pCorrect = nan(1,1);
                    ta(iarea).perm(isubpop).dir(idir).mean_Error = nan(1,1);
                    ta(iarea).perm(isubpop).dir(idir).mean_bias = nan(1,1);
                    end
                end
            else % if enough trials, do decoding

                for iarea = 1:numel(areas)

                    tempAreaUnits = goodUnits(strcmp([goodUnits.ecephys_structure_acronym],areas(iarea)));
                    nUnits = numel(tempAreaUnits);
                    nSubPops = floor(nUnits/reqUnits); % number of subpops we can gert

                    if nSubPops<1 % if not enough units to do any decoding, set to nan
                        for idir2=1:4
                        ta(iarea).perm(1).dir(idir2).mean_pCorrect = nan(1,1);
                        ta(iarea).perm(1).dir(idir2).mean_Error = nan(1,1);
                        ta(iarea).perm(1).dir(idir2).mean_bias = nan(1,1);
                        end
                    else % if enough units, do decoding

                        tempAreaUnits = tempAreaUnits(randperm(nUnits)); % shuffle order of units for the sake of it

                        for isubpop = 1:nSubPops
                            tempUnits = tempAreaUnits(isubpop*reqUnits-(reqUnits-1):(isubpop*reqUnits)); % units for this subpop

                            % get the trials to be used for decoding (downsampling to nTrials)
                            for icond = 1:nConds
                                nUseableTrials = size(tempUnits(1).(scName){idir,icond},1);
                                randSampTrials = randsample(nUseableTrials,reqTrials);
                                cond{icond}.trials = randSampTrials;
                            end


                            % generate input struct for decoding function
                            units2Decode = []; units2Decode(reqUnits).spikeCell = [];

                            for iunit = 1:numel(tempUnits)
                                for icond = 1:nConds
                                    units2Decode(iunit).spikeCell{icond} =...
                                        tempUnits(iunit).(scName){idir,icond}...
                                        (cond{icond}.trials);
                                end
                            end

                            % do the LOO-CV decoding
                            tuningAxes = {[1:7]};
                            ta(iarea).perm(isubpop).dir(idir).cond = PID_DecodeFromSpikeCellArray(units2Decode, tuningAxes, 0);

                            % get pCorrect and bias of predictions for each speed, for
                            % this subpop
                            for icond = 1:nConds
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).pCorrect = prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds==icond);
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).meanError = mean(abs(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds-icond));
                                ta(iarea).perm(isubpop).dir(idir).cond(icond).bias = prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds>icond) - prop(ta(iarea).perm(isubpop).dir(idir).cond(icond).preds<icond);
                            end

                            ta(iarea).perm(isubpop).dir(idir).mean_pCorrect = mean([ta(iarea).perm(isubpop).dir(idir).cond.pCorrect]);
                            ta(iarea).perm(isubpop).dir(idir).mean_Error = mean([ta(iarea).perm(isubpop).dir(idir).cond.meanError]);
                            ta(iarea).perm(isubpop).dir(idir).mean_bias = mean([ta(iarea).perm(isubpop).dir(idir).cond.bias]);
                        end

                    end

                end

            end

        end

    else % set all nans if no trials
        for iarea = 1:8
            for idir = 1:4
                ta(iarea).perm(1).dir(idir).mean_pCorrect = nan(1,1);
                ta(iarea).perm(1).dir(idir).mean_Error = nan(1,1);
                ta(iarea).perm(1).dir(idir).mean_bias = nan(1,1);
            end
        end

    end

    run.ta=ta;
    clear ta;



    %% save data

    save(outputFileName, 'stat','run', 'sessionInfo')

end

end
