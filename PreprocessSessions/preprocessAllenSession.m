function preprocessAllenSession(sessionID, inputFile, outputFile)

load(inputFile)

%%
%%%% process units %%%%

%% generate units struct with spike times and generic info
units = genUnitStructFromAllenData(unitInformation,unitIndex,spikeTimes);

%% get layer of cortical units based on CCF co-ordinate
% this uses AllenCCF github repo from Cortexlab to load the CCF
coordsArray = [vertcat(units.anterior_posterior_ccf_coordinate),...
    vertcat(units.dorsal_ventral_ccf_coordinate),...
    vertcat(units.left_right_ccf_coordinate)];

[~, ~, ~, layer] = getAllenAcronymFromCCF(coordsArray);

%
for iunit=1:numel(units)
    units(iunit).genotype = {sessionInfo.genotype};
    units(iunit).sessionID = sessionID;
    units(iunit).layer = layer(iunit);
end

%% check if units are optotagged as PV, SST or VIP interneurons
units(iunit).optotag = deal(nan);
if strcmp(sessionInfo.genotype, 'wt/wt')
elseif strcmp(sessionInfo.genotype, 'Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt') %PV optotag
    load('PVTags.mat')
    idx = find(ismember([units.ID],PVtags));
    for iunit = 1:numel(idx)
        units(idx(iunit)).optotag = 'PV';
    end
elseif strcmp(sessionInfo.genotype, 'Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt') % SST optotag
    load('SstTags.mat')
    idx = find(ismember([units.ID],SSTtags));
    for iunit = 1:numel(idx)
        units(idx(iunit)).optotag = 'SST';
    end
elseif strcmp(sessionInfo.genotype, 'Vip-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt') % VIP optotag
    load('VIPtags.mat')
    idx = find(ismember([units.ID],VIPtags));
    for iunit = 1:numel(idx)
        units(idx(iunit)).optotag = 'VIP';
    end
end


%%
%%%% process trials %%%%

% time intervals around the start and stop time of a stimulus we are
% interested in.
pretime = 0.2;
posttime = 0.2;

%% generate trial struct

trials = genTrialStructFromAllenData(stim_table, false);

%% tag trials with invalid times (probe errors)
[trials.containsInvalidTime] = deal(false);

if isfield(sessionInfo,'invalid_times')
for irow = 1:size(sessionInfo.invalid_times,1)
    interval(1) = sessionInfo.invalid_times{irow,1};
    interval(2) = sessionInfo.invalid_times{irow,2};
    
    % find trials that start or finish within the invalid time.
    idx = find( ([trials.start_time]-pretime<interval(1) & [trials.stop_time]+posttime>interval(1)) |...
        ([trials.stop_time]+posttime>interval(1) & [trials.start_time]-pretime<interval(2)));
    try
    [trials(idx).containsInvalidTime] = deal(true);
    catch
        debug = 1;
    end

end
end

%% add run and pupil info to trials struct

% resample and smooth run trace with 50ms moving average.
samplingFreq = 100;
smthWindow = 0.175;

smthWindowBin = smthWindow/(1/samplingFreq);
sampleSpacing = 1/samplingFreq;

xq = min(runInfo.runMidTimes):sampleSpacing:max(runInfo.runMidTimes); % new interp time
runInfo.runTimeInterp = xq;
runInfo.runSpeedInterp = interp1(runInfo.runMidTimes, runInfo.runSpeed, xq); 
runInfo.smthRunSpeedInterp = smoothdata(runInfo.runSpeedInterp,'gaussian',smthWindowBin);


trials = getRunAndPupilInfo(trials, runInfo, pupilInfo, pretime, posttime);


%% save data
save(outputFile, 'units', 'trials', 'runInfo', 'pupilInfo', 'sessionInfo')


end