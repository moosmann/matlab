classdef GenericProfileWizardManager
    %GENERICPROFILEWIZARDMANAGER - Manager for the JavaScript Generic Profile Wizard
    
    %   Copyright 2017 The MathWorks, Inc.
    methods (Hidden, Static)
        
        function installedSchedulers = getInstalledSchedulers()
            % Retrieve information about where support packages are installed
            installationFolder = fullfile(matlabshared.supportpkg.getSupportPackageRoot(), 'parallel', '{0}');
            integrationScriptsLocation = fullfile(installationFolder, '{1}');
            matlabVersion = ['R' version('-release')];
            
            installedSchedulers = struct(...
                'installationFolder', installationFolder, ...
                'integrationScriptsLocation', integrationScriptsLocation, ...
                'matlabVersion', matlabVersion);
            
            % Determine which scheduler support packages are installed
            allInstalledSupportPackages = matlabshared.supportpkg.getInstalled();
            for supportPkg = allInstalledSupportPackages
                try
                    sched = parallel.internal.types.SupportPackageScheduler.fromName(supportPkg.Name);
                    installedSchedulers.(sched.FolderName) = supportPkg;
                catch
                    % Not a PCT support package for integration with a 3rd party scheduler
                end
            end
        end
        
        function profileName = createProfile(model)
            c = parallel.cluster.Generic();
            
            % Scheduler and submission mode configuration
            c.OperatingSystem = model.operatingSystem;
            c.IntegrationScriptsLocation = model.integrationScriptsLocation;
            submissionMode = parallel.internal.types.SubmissionMode.fromName(model.submissionMode);
            
            if submissionMode.RequiresClusterHost
                c.AdditionalProperties.ClusterHost = model.clusterHost;
            end
            
            if submissionMode.RequiresRemoteJobStorageLocation
                c.AdditionalProperties.RemoteJobStorageLocation = model.remoteJobStorageLocation;
                c.AdditionalProperties.UseUniqueSubfolders = model.useUniqueSubfolders;
            end
            
            c.HasSharedFilesystem = model.hasSharedFileSystem;
            
            % Other configuration
            c.NumWorkers = model.numWorkers;
            c.ClusterMatlabRoot = model.clusterMatlabRoot;
            c.NumThreads = model.numThreads;
            c.RequiresMathWorksHostedLicensing = strcmp(model.licenseType, 'mhlm');
            if c.RequiresMathWorksHostedLicensing
                c.hSetPropertyNoCheck('LicenseNumber', model.licenseNumber);
            end
            
            % Save profile with the given name and (optional) description
            profileName = model.profileName;
            clusterDescription = model.clusterDescription;
            
            if isempty(clusterDescription)
                c.saveAsProfile(profileName);
            else
                c.saveAsProfile(profileName, 'Description', clusterDescription);
            end
        end
        
        function openClusterProfileManager()
            com.mathworks.toolbox.distcomp.ui.desk.ProfileDesktopMenuProvider.showProfileManagerUI();
        end
        
        function openDocPage(pathToMapFileFromDocroot, anchor)
            helpview(fullfile(docroot, pathToMapFileFromDocroot), anchor);
        end
    end
end
