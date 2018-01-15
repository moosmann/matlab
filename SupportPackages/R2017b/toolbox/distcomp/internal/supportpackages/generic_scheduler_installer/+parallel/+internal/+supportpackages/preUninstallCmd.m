function preUninstallCmd(schedulerName, relativeScriptPath)
% PREUNINSTALLCMD - internal: performs generic scheduler uninstallation workflow

% Copyright 2015-2017 The MathWorks, Inc.

installationFolder = matlabshared.supportpkg.getSupportPackageRoot();
parallel.internal.supportpackages.registerMessageCatalog(installationFolder);

iPreUninstallBackupWorkflow(installationFolder, installationFolder, schedulerName, relativeScriptPath);
end

function iPreUninstallBackupWorkflow(installationFolder, installRoot, schedulerName, relativeScriptPath)
iLogMessage(sprintf('preUninstallBackupWorkflow(%s, %s, %s)', installationFolder, schedulerName, relativeScriptPath));
sourceLocation = fullfile(installationFolder, relativeScriptPath);
defaultBackupLocation = pwd;
defaultName = iGetDefaultBackupFilename(schedulerName);
destinationLocation = fullfile(defaultBackupLocation, defaultName);
state = 'askingIfBackupRequired';
while true
    iLogMessage(sprintf('preUninstallBackupWorkflow state: %s', state));
    switch state
        case 'askingIfBackupRequired'
            [cancelled, backupRequired] = iAskUserIfBackupRequired(sourceLocation, schedulerName);
            if cancelled
                % We cannot abort uninstallation once we're in preUninstallCmd due to SPI constraints.
                % Therefore we ask user again until they provide an answer to the backup question.
                state = 'askingIfBackupRequired';
            elseif backupRequired
                state = 'choosingLocation';
            else
                state = 'commenceUninstall';
            end
        case 'choosingLocation'
            [cancelled, destinationLocation] = iAskUserForBackupLocation(defaultBackupLocation, defaultName);
            if cancelled
                state = 'askingIfBackupRequired';
            elseif ~iBackupLocationValid(installRoot, destinationLocation)
                state = 'badBackupPath';
            else
                state = 'performingBackup';
            end
        case 'badBackupPath'
            iBadBackupLocationDialog();
            state = 'choosingLocation';
        case 'performingBackup'
            try
                iPerformBackup(destinationLocation, sourceLocation);
            catch exception
                iBackupFailureDialog(destinationLocation, sourceLocation);
                rethrow(exception);
            end
            state = 'backupSuccess';
        case 'backupSuccess'
            iBackupSuccessDialog(destinationLocation);
            state = 'commenceUninstall';
        case 'commenceUninstall'
            break;
        otherwise
            error('Undefined state');
    end
end
end

function [cancelled, backupRequired] = iAskUserIfBackupRequired(sourceLocation, schedulerName)
iLogMessage(sprintf('askUserIfBackupRequired(%s, %s)', sourceLocation, schedulerName));

backupQueryMessage = getString(message('parallel_supportpackages:generic_scheduler:UninstallBackupQueryMessage', schedulerName, sourceLocation));
backupQueryTitle = getString(message('parallel_supportpackages:generic_scheduler:UninstallBackupQueryTitle'));

maximumRetries = 10;
failureCount = 0;
finished = false;
while ~finished && failureCount < maximumRetries
    try
        result = questdlg(backupQueryMessage , backupQueryTitle, 'Yes', 'No', 'Yes');
        % result is empty if user closes questdlg window
        finished = true;
    catch exception
        % If uncaught, any exceptions from questdlg would cause the uninstallation to
        % commence without the user choosing to create or skip a backup.
        iLogMessage(sprintf('askUserIfBackupRequired exception: %s', exception.identifier));
        if strcmp(exception.identifier, 'MATLAB:hg:DeletedObject')
            % If questdlg is closed by the user before it has fully constructed,
            % an MATLAB:hg:DeletedObject exception is thrown. Retry again.
            failureCount = failureCount + 1;
        else
            % This is an unknown exception - we stop retrying without finishing.
            break;
        end 
    end
end
if finished
    cancelled = isempty(result);
    backupRequired = strcmp('Yes', result);
else
    % We were unable to finish successfully,
    % so we take the most conservative action possible, and assume the user wants a backup.
    cancelled = false;
    backupRequired = true;
end
end

function name = iGetSafeSchedulerFilename(schedulerName)
lowerCaseSchedulerName = lower(schedulerName);
name = strrep(lowerCaseSchedulerName, ' ', '');
end

function defaultName = iGetDefaultBackupFilename(schedulerName)
iLogMessage(sprintf('getDefaultBackupFilename(%s)', schedulerName));

safeSchedulerName = iGetSafeSchedulerFilename(schedulerName);
defaultName = strcat(safeSchedulerName, '_backup.zip');
end

function [cancelled, fullLocation] = iAskUserForBackupLocation(defaultLocation, defaultName)
iLogMessage(sprintf('askUserForBackupLocation(%s, %s)', defaultLocation, defaultName));

backupLocationQueryTitle = getString(message('parallel_supportpackages:generic_scheduler:BackupLocationQueryTitle'));
fullLocation = fullfile(defaultLocation, '*.zip');
[filename, location] = uiputfile(fullLocation, backupLocationQueryTitle, defaultName);
% uiputfile sets filename, location to 0 if user presses cancel
cancelled = isequal(filename,0) || isequal(location,0);
fullLocation = fullfile(location, filename);
end

function locationValid = iBackupLocationValid(sourceLocation, destinationLocation)
iLogMessage(sprintf('backupLocationValid(%s, %s)', sourceLocation, destinationLocation));

% Do not permit user to choose installation location as backup location
locationValid = ~strncmp(sourceLocation, destinationLocation, length(sourceLocation));
end

function iBadBackupLocationDialog()
iLogMessage('badBackupLocationDialog()');

warningMessage = getString(message('parallel_supportpackages:generic_scheduler:BadBackupLocationMessage'));
title = getString(message('parallel_supportpackages:generic_scheduler:BadBackupLocationTitle'));
dlg = warndlg(warningMessage, title, 'modal');
uiwait(dlg);
end

function iPerformBackup(destination, source)
iLogMessage(sprintf('performBackup(%s,%s)', destination, source));

zip(destination, source);
end

function iBackupSuccessDialog(destination)
iLogMessage(sprintf('backupSuccessDialog(%s)', destination));

successMessage = getString(message('parallel_supportpackages:generic_scheduler:BackupSuccessMessage', destination));
title = getString(message('parallel_supportpackages:generic_scheduler:BackupSuccessTitle'));
dlg = helpdlg(successMessage, title);
uiwait(dlg);
end

function iBackupFailureDialog(destination, source)
iLogMessage(sprintf('backupFailureDialog(%s, %s)', destination, source));

errorMessage = getString(message('parallel_supportpackages:generic_scheduler:BackupFailureMessage', destination, source));
title = getString(message('parallel_supportpackages:generic_scheduler:BackupFailureTitle'));
dlg = errordlg(errorMessage, title);
uiwait(dlg);
end

function iLogMessage(message)
parallel.internal.supportpackages.logMessage(message);
end
