%% create input and output args for python script

% folder = 'D:\AllenSDK\FunctionalConnectivity\';
% 
% % get data folders
% nwbfolders = dir(fullfile(folder, 'session_*'));
% dirFlags = [nwbfolders.isdir];
% nwbfolders = nwbfolders(dirFlags);
% 
% prefix = 'session_';  % The prefix to remove
% regex = {['^' prefix]}; % the regular expressions for prefix & suffix
% replacements = {''}; % Replacement strings
% 
% failedFiles = [];
load('FCsessions.mat')

%
%%
for k = 1 : numel(fc_sessionID)
  pyargs(k).sessionID = string(fc_sessionID(k));
  %pyargs(k).outputFile = ['C:\Users\edward.horrocks\Documents\GitHub\AllenDataAnalysis\Data\rawMatFiles\session_' num2str(pyargs(k).sessionID) '.mat'];
  pyargs(k).outputFile = ['D:\AllenSDK\FunctionalConnectivity22\session_' num2str(pyargs(k).sessionID) '.mat'];
  %pyargs(k).sysCommand = ['C:\Users\edward.horrocks\PycharmProjects\allen-data\venv\Scripts\python.exe C:\Users\edward.horrocks\PycharmProjects\allen-data\venv\getCSDplots.py -i ' char(pyargs(k).sessionID)];
  %pyargs(k).sysCommand = ['C:\Users\edward.horrocks\PycharmProjects\allen-data\venv\Scripts\python.exe C:\Users\edward.horrocks\PycharmProjects\allen-data\venv\allenSession2matfile.py -i ' char(pyargs(k).sessionID) ' -o ' char(pyargs(k).outputFile)];
  pyargs(k).sysCommand = ['C:\Users\edward.horrocks\PycharmProjects\allen2\venv\Scripts\python.exe C:\Users\edward.horrocks\PycharmProjects\allen2\venv\allenSession2matfile.py -i ' char(pyargs(k).sessionID) ' -o ' char(pyargs(k).outputFile)];

end

%%
tic

for k = 1:numel(pyargs)
    k
    try
    system(pyargs(k).sysCommand)
    catch
        fprintf(['failed to process session: k = ' num2str(k) '\n'])
        failedFiles = [failedFiles, k];
    end
end
toc

%%
    
    
