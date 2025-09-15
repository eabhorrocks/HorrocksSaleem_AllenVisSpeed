%% Batch processing

%% get session filenames

dataDir = 'D:\AllenMatFiles\FC';
% dataFiles = dir(fullfile(dataDir, 'pre_session*'));

for ifile = 1:numel(dataFiles)
    dataFiles(ifile).sessionID = str2double(regexp(dataFiles(ifile).name, '\d*', 'match'));
    dataFiles(ifile).filename = fullfile(dataFiles(ifile).folder, dataFiles(ifile).name);
end

%% get session filenames

savePrefix = 'FC_processed';

tic
parfor isession = 1:numel(dataFiles)
    
isession
allen_analyseSesssion_DotMotionSpeedTuning(dataFiles(isession).filename, fullfile(dataFiles(isession).folder, [savePrefix,'_', num2str(dataFiles(isession).sessionID), '.mat']))
allen_analyseSesssion_TuningDGs(dataFiles(isession).filename, fullfile(dataFiles(isession).folder, [savePrefix,'_', num2str(dataFiles(isession).sessionID), '.mat']))
allen_analyseSesssion_PIDSpeedTuning(dataFiles(isession).filename, fullfile(dataFiles(isession).folder, [savePrefix,'_', num2str(dataFiles(isession).sessionID), '.mat']))

end
toc
