function postInstallCmd(scheduler)
%POSTINSTALLCMD - internal: displays generic scheduler post installation dialog

% Copyright 2015-2017 The MathWorks, Inc.
validateattributes(scheduler, {'parallel.internal.types.SupportPackageScheduler'}, {'scalar'});

installationFolder = matlabshared.supportpkg.getSupportPackageRoot();
parallel.internal.supportpackages.registerMessageCatalog(installationFolder);

displayName = scheduler.DisplayName;
folderName = scheduler.FolderName;
scriptLocation = fullfile(installationFolder, 'parallel', folderName);

dialogMessage = getString(message('parallel_supportpackages:generic_scheduler:InstallationCompleteMessage', ...
    displayName, scriptLocation));
dialogTitle = getString(message('parallel_supportpackages:generic_scheduler:InstallationCompleteTitle'));

dialog = msgbox(dialogMessage, dialogTitle, 'modal');
uiwait(dialog);

parallel.cluster.generic.runProfileWizard(scheduler);
end
