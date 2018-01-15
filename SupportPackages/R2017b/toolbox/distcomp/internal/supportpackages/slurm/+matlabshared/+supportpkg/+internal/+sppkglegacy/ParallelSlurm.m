%Support package registry plugin for SLURM
%   This class is used by support package infrastructure to find the legacy
%   support_package_registry.xml for the BaseCode specified
%
%   Type:
%   >> help matlabshared.supportpkg.internal.sppkglegacy.SupportPackageRegistryPluginBase
%   for more information.

%   Copyright 2016 The MathWorks, Inc.

classdef ParallelSlurm < matlabshared.supportpkg.internal.sppkglegacy.SupportPackageRegistryPluginBase
    
    properties(Constant)
        BaseCode = 'PCT_SLURM';
    end
end
