%% load data

toAnalyse = 'BO'; % 'FC' or 'BO'

if strcmp(toAnalyse, 'FC')
    dataDir = 'D:\AllenMatFiles\FC';
    dataFiles = dir(fullfile(dataDir, 'FC_processsed*')); %
elseif strcmp(toAnalyse, 'BO')
    dataDir = 'D:\AllenMatFiles\BO';
    dataFiles = dir(fullfile(dataDir, 'BO_processed*')); %
end


tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'units');
    session(ifile).units = units;
end
toc
%% generate allUnits struct

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};

allUnits = [session.units];

nUnits = numel(allUnits);

goodUnits = allUnits([allUnits.isi_violations]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.waveform_amplitude]>=50 &...
    ismember([allUnits.ecephys_structure_acronym],areas));
nGoodUnits = numel(goodUnits);