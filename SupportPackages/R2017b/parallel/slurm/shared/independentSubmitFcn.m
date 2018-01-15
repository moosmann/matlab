function independentSubmitFcn(cluster, job, environmentProperties)
%INDEPENDENTSUBMITFCN Submit a MATLAB job to a Slurm cluster
%
% Set your cluster's IntegrationScriptsLocation to the parent folder of this
% function to run it when you submit an independent job.
%
% See also parallel.cluster.generic.independentDecodeFcn.
%

% Copyright 2010-2016 The MathWorks, Inc.

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;
if ~isa(cluster, 'parallel.Cluster')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The function %s is for use with clusters created using the parcluster command.', currFilename)
end

decodeFunction = 'parallel.cluster.generic.independentDecodeFcn';

if ~cluster.HasSharedFilesystem
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end


if ~strcmpi(cluster.OperatingSystem, 'unix')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s only supports clusters with unix OS.', currFilename)
end

% The job specific environment variables
% Remove leading and trailing whitespace from the MATLAB arguments
matlabArguments = strtrim(environmentProperties.MatlabArguments);
variables = {'MDCE_DECODE_FUNCTION', decodeFunction; ...
    'MDCE_STORAGE_CONSTRUCTOR', environmentProperties.StorageConstructor; ...
    'MDCE_JOB_LOCATION', environmentProperties.JobLocation; ...
    'MDCE_MATLAB_EXE', environmentProperties.MatlabExecutable; ...
    'MDCE_MATLAB_ARGS', matlabArguments; ...
    'MDCE_DEBUG', 'true'; ...
    'MLM_WEB_LICENSE', environmentProperties.UseMathworksHostedLicensing; ...
    'MLM_WEB_USER_CRED', environmentProperties.UserToken; ...
    'MLM_WEB_ID', environmentProperties.LicenseWebID; ...
    'MDCE_LICENSE_NUMBER', environmentProperties.LicenseNumber; ...
    'MDCE_STORAGE_LOCATION', environmentProperties.StorageLocation};
% Set each environment variable to newValue if currentValue differs.
% We must do this particularly when newValue is an empty value,
% to be sure that we clear out old values from the environment.
for ii = 1:size(variables, 1)
    variableName = variables{ii,1};
    currentValue = getenv(variableName);
    newValue = variables{ii,2};
    if ~strcmp(currentValue, newValue)
        setenv(variableName, newValue);
    end
end

% Deduce the correct quote to use based on the OS of the current machine
if ispc
    quote = '"';
else
    quote = '''';
end


% The script name is independentJobWrapper.sh
scriptName = 'independentJobWrapper.sh';
% The wrapper script is in the same directory as this file
dirpart = fileparts(mfilename('fullpath'));
quotedScriptName = sprintf('%s%s%s', quote, fullfile(dirpart, scriptName), quote);

% Get the tasks for use in the loop
tasks = job.Tasks;
numberOfTasks = environmentProperties.NumberOfTasks;
jobIDs = cell(numberOfTasks, 1);
% Loop over every task we have been asked to submit
for ii = 1:numberOfTasks
    taskLocation = environmentProperties.TaskLocations{ii};
    % Set the environment variable that defines the location of this task
    setenv('MDCE_TASK_LOCATION', taskLocation);
    
    % Choose a file for the output. Please note that currently, JobStorageLocation refers
    % to a directory on disk, but this may change in the future.
    logFile = cluster.getLogLocation(tasks(ii));
    quotedLogFile = sprintf('%s%s%s', quote, logFile, quote);
    
    % Submit one task at a time
    jobName = sprintf('Job%d.%d', job.ID, tasks(ii).ID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CUSTOMIZATION MAY BE REQUIRED %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    additionalSubmitArgs = '--ntasks=1';
    commonSubmitArgs = getCommonSubmitArgs(cluster);
    if ~isempty(commonSubmitArgs) && ischar(commonSubmitArgs)
        additionalSubmitArgs = strtrim([additionalSubmitArgs, ' ', commonSubmitArgs]);
    end
    dctSchedulerMessage(5, '%s: Generating command for task %i', currFilename, ii);
    commandToRun = getSubmitString(jobName, quotedLogFile, quotedScriptName, ...
        additionalSubmitArgs);
    
    % Now ask the cluster to run the submission command
    dctSchedulerMessage(4, '%s: Submitting job using command:\n\t%s', currFilename, commandToRun);
    try
        % Make the shelled out call to run the command.
        [cmdFailed, cmdOut] = system(commandToRun);
    catch err
        cmdFailed = true;
        cmdOut = err.message;
    end
    if cmdFailed
        error('parallelexamples:GenericSLURM:SubmissionFailed', ...
            'Submit failed with the following message:\n%s', cmdOut);
    end
    
    dctSchedulerMessage(1, '%s: Job output will be written to: %s\nSubmission output: %s\n', currFilename, logFile, cmdOut);
    jobIDs{ii} = extractJobId(cmdOut);
    
    if isempty(jobIDs{ii})
        warning('parallelexamples:GenericSLURM:FailedToParseSubmissionOutput', ...
            'Failed to parse the job identifier from the submission output: "%s"', ...
            cmdOut);
    end
end

% set the job ID on the job cluster data
cluster.setJobClusterData(job, struct('ClusterJobIDs', {jobIDs}));
