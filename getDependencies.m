%% Get external, non-native matlab dependencies
% Finds, saves, and copies only EXTERNAL non-native function dependencies.

clear; clc;

% --- User Configuration ---
outputCsvFile = 'dependency_report.csv';
dependencyFolderName = 'imported_dependencies'; % Folder to copy files into

%% PART 1: Find Dependencies (Manual Scan Method)
disp('--- PART 1: Finding all non-native dependencies ---');
allFiles = dir('**/*.m');
allFiles(startsWith({allFiles.folder}, fullfile(pwd, '.git'))) = [];
disp(['Found ' num2str(length(allFiles)) ' files to scan.']);

allPotentialFunctions = {};
for i = 1:length(allFiles)
    filePath = fullfile(allFiles(i).folder, allFiles(i).name);
    try
        content = fileread(filePath);
        [matchStart, matchEnd] = regexp(content, '[a-zA-Z]\w*');
        for j = 1:length(matchStart)
            word = content(matchStart(j):matchEnd(j));
            nextCharIndex = matchEnd(j) + 1;
            while nextCharIndex <= length(content) && isspace(content(nextCharIndex))
                nextCharIndex = nextCharIndex + 1;
            end
            if nextCharIndex <= length(content) && content(nextCharIndex) == '('
                allPotentialFunctions = [allPotentialFunctions; {word}];
            end
        end
    catch ME
        fprintf('Could not read or parse file: %s (%s)\n', filePath, ME.message);
    end
end

uniqueNames = unique(allPotentialFunctions);
disp(['Found ' num2str(length(uniqueNames)) ' unique potential function names.']);

matlabRoot = matlabroot;
nonNativeFunctions = {};
matlabKeywords = {'if', 'for', 'while', 'switch', 'case', 'else', 'elseif', 'end', 'try', 'catch', 'function', 'return', 'parfor', 'classdef', 'properties', 'methods', 'events'};

for i = 1:length(uniqueNames)
    funcName = uniqueNames{i};
    if ismember(funcName, matlabKeywords)
        continue;
    end
    try
        location = which(funcName);
        
        % --- FIX #1: More robust check for built-in and native functions ---
        isBuiltIn = startsWith(location, 'built-in', 'IgnoreCase', true);
        isNative = startsWith(location, matlabRoot);
        
        if ~isempty(location) && ~isBuiltIn && ~isNative && ~strcmpi(location, 'variable')
             nonNativeFunctions = [nonNativeFunctions; {funcName, location}];
        end
        % --------------------------------------------------------------------
    catch
        continue;
    end
end
disp('Dependency analysis complete.');
fprintf('\n');

% --- Central Checkpoint ---
if isempty(nonNativeFunctions)
    disp('No non-native dependencies found. Nothing to save or copy.');
    disp('--- All tasks finished. ---');
    return;
end

%% PART 2: Save Results to a Table
disp('--- PART 2: Saving results to a table ---');
resultsTable = cell2table(nonNativeFunctions, 'VariableNames', {'FunctionName', 'SourcePath'});
writetable(resultsTable, outputCsvFile);
disp(['Results saved to "' outputCsvFile '"']);
fprintf('\n');

%% PART 3: Copy External Dependency Files
disp('--- PART 3: Copying EXTERNAL dependency files ---');
if ~exist(dependencyFolderName, 'dir')
    mkdir(dependencyFolderName);
end

disp(['Copying files to "' dependencyFolderName '" while preserving structure...']);
errorCount = 0;
copyCount = 0;
repoDir = pwd;

for i = 1:height(resultsTable)
    try
        sourcePath = resultsTable.SourcePath{i};
        [sourceDir, fileName, fileExt] = fileparts(sourcePath);
        
        % --- FIX #2: Add a failsafe to ensure the source is a valid file ---
        if ~isfile(sourcePath)
            fprintf('Skipping (not a valid file): %s\n', sourcePath);
            continue;
        end
        % --------------------------------------------------------------------
        
        if startsWith(sourcePath, repoDir)
            fprintf('Skipping (already in repo): %s\n', [fileName, fileExt]);
            continue;
        end

        pathParts = split(sourceDir, filesep);
        parentFolder = pathParts{end};
        
        sanitizedParentFolder = regexprep(parentFolder, '[\\/:"*?<>|]', '_');
        sanitizedFileName = regexprep(fileName, '[\\/:"*?<>|]', '_');
        
        destinationDir = fullfile(pwd, dependencyFolderName, sanitizedParentFolder);
        destinationPath = fullfile(destinationDir, [sanitizedFileName, fileExt]);
        
        if ~exist(destinationDir, 'dir')
            mkdir(destinationDir);
        end
        
        copyfile(sourcePath, destinationPath);
        fprintf('Copied: %s\n', [fileName, fileExt]);
        copyCount = copyCount + 1;
        
    catch ME
        fprintf('ERROR copying "%s": %s\n', [fileName, fileExt], ME.message);
        errorCount = errorCount + 1;
    end
end

fprintf('\nFile copying complete. Copied %d external files.\n', copyCount);
if errorCount > 0
    fprintf('%d errors occurred during the file copy process.\n', errorCount);
end

disp('--- All tasks finished. ---');