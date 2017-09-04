% Batch Processing Tool
% BatchProcessing.m
%
%   This code allows the user to select multiple patient folders to pass to
%   main.m for processing. It will attempt to run main.m on all selected
%   patients, but skip patients that throw errors. A log will be saved. If
%   no second argument is passed to main.m, the results will be saved in
%   each patient's folder, but if a results folder is passed, all results
%   will be saved in a single folder.
%
%   Inputs: 
%   - User selects patient folders. Each patient folder should have a
%   subfolder, "NUNDA" that contains a structure with mask and other data,
%   called "NUNDAin.mat"
%
%   Constants:
%
%   Outputs:
%   - log: a diary file will be saved in the directory noted at the end of
%   execution.
%   - NUNDAout.mat: a structure containing the results of the calculations
%
%   Example usage:
%   - script, does not accept input arguments
%
%   Required functions:
%   - main.m and all dependencies
%   - uipickfiles: written by Douglas M. Schwarz, from Mathworks file
%   exchange
%
%   References:
%
%   Known bugs/shortcomings:
%   - Add a switch to change if all results are compiled in one folder or
%   in patient folders
%
% Written by Mike Scott, August 2017 (Northwestern University)
% michael.scott1@northwestern.edu

addpath(genpath([pwd() filesep() 'functions']));

% Get the folders that should be processed
[folders] = uipickfiles('Prompt','Select all patient folders to process','FilterSpec','C:\Users\Mike\Desktop\Working');

% Start a log file
clc; close all;

% Make a folder to store the results in the parent directory
[folder_name,~,~] = fileparts(folders{1});
formatOut = 'yyyymmdd_HHMMSS';
results_folder = [folder_name filesep() datestr(now,formatOut) '_results'];
mkdir(results_folder)
diary([results_folder filesep() 'log.txt']);
fprintf('==================================================\n');
fprintf('===           Batch Processing Tool            ===\n');
fprintf('==================================================\n\n');
overalltime = tic;
fprintf('%i folders chosen, starting execution...\n\n',length(folders));

% Iterate through all the chosen folders
counter = 0;
for ii = 1:length(folders)
    fprintf('Beginning patient %i/%i...\n\n',ii,length(folders))
    % Get the patient folder name to use in the results folder
    [~,patientID,~] = fileparts(folders{ii});
    try
        cardiac4dflow(folders{ii},[results_folder filesep() patientID]);
        close all;
    catch
        warning('Analysis had an error in folder:\n   %s',folders{ii});
        counter = counter + 1;
    end
end

fprintf('==================================================\n');
fprintf('==================================================\n');
fprintf('==================================================\n\n');
fprintf('Code execution complete in %i seconds.\n',round(toc(overalltime)));
fprintf('   Code ran with errors in %i subjects.\n',counter)
fprintf('   Execution log saved to:\n     %s\n',[results_folder filesep() 'log.txt'])
fprintf('================================================:)\n\n');
diary off



