%% load sessions

dataDir = 'D:\AllenMatFiles\FC'; 
%dataDir = 'D:\AllenMatFiles\BO';


dataFiles = dir(fullfile(dataDir, 'session*.mat'));


%%
savePrefix = 'pre_session_';

for ifile =1:numel(dataFiles)
    ifile
    sessionID = str2double(regexp(dataFiles(todo(ifile)).name, '\d*', 'match'));
    inputFileName = fullfile(dataDir, dataFiles(todo(ifile)).name);
    outputFileName = fullfile(dataDir, [savePrefix, num2str(sessionID), '.mat']);

    preprocessAllenSession(sessionID, inputFileName, outputFileName)
end