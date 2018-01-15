% SupportPackageScheduler - Enumeration of 3rd party schedulers for which PCT provides integration support packages

%   Copyright 2017 The MathWorks, Inc.

classdef SupportPackageScheduler < parallel.internal.types.NamedEnumeration
    enumeration
        GridEngine('Grid Engine', 'gridengine')
        LSF('IBM Platform LSF', 'lsf')
        PBS('PBS', 'pbs')
        Slurm('Slurm', 'slurm')
    end
    
    properties (Constant, GetAccess = private)
        SupportPackageNamePrefix = 'Parallel Computing Toolbox integration for MATLAB Distributed Computing Server with ';
    end
    
    properties (SetAccess = immutable)
        DisplayName
        FolderName
        SupportPackageName
    end
    
    methods
        function obj = SupportPackageScheduler(displayName, folderName)
            supportPackageName = [parallel.internal.types.SupportPackageScheduler.SupportPackageNamePrefix, displayName];
            obj@parallel.internal.types.NamedEnumeration({displayName, folderName, supportPackageName});
            obj.DisplayName = displayName;
            obj.FolderName = folderName;
            obj.SupportPackageName = supportPackageName;
        end
    end
    
    methods (Static)
        function enumVal = fromName(name)
            enumVal = parallel.internal.types.NamedEnumeration.getEnumFromName(...
                'parallel.internal.types.SupportPackageScheduler', name);
        end
    end
end