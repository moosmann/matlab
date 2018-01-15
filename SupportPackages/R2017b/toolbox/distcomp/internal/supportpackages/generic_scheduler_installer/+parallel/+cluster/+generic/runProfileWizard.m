function runProfileWizard(scheduler)
%RUNPROFILEWIZARD Runs the wizard for creating a new Generic profile

% Copyright 2017 The MathWorks, Inc.
    schedulerSpecified = nargin > 0;
    if schedulerSpecified
        validateattributes(scheduler, {'parallel.internal.types.SupportPackageScheduler'}, {'scalar'});
    end

    % Ensure message catalog is registered
    installationFolder = matlabshared.supportpkg.getSupportPackageRoot();
    parallel.internal.supportpackages.registerMessageCatalog(installationFolder);

    % Create the dialog window
    connector.ensureServiceOn();
    addOnName = 'pctgenericprofilewizard';
    connector.addWebAddOnsPath(addOnName, [fileparts(mfilename('fullpath')) ...
        '/../../../../generic_scheduler_installer_web']);
    url = connector.getUrl(sprintf('addons/%s/index.html', addOnName));
    if schedulerSpecified
        schedulerKey = scheduler.FolderName;
        url = [url '&scheduler=' schedulerKey];
    end
    wizardWindow = matlab.internal.webwindow(url, 4440);
    if schedulerSpecified
        displayName = scheduler.DisplayName;
        wizardWindow.Title = getString(message('parallel_supportpackages:generic_scheduler:WizardTitlePostInstall', displayName));
    else
        wizardWindow.Title = getString(message('parallel_supportpackages:generic_scheduler:WizardTitle'));
    end

    % Adjust the size and position of the wizard
    screenSize = get(0, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);
    width = 640;
    height = 440;
    wizardWindow.Position = [(screenWidth - width) / 2, (screenHeight - height) / 2, width, height];
    wizardWindow.setMinSize([width height]);

    % Show the wizard
    wizardWindow.show();
end