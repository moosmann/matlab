function postInstallCmd()
%POSTINSTALLCMD - internal: post install hook for Slurm

% Copyright 2015-2017 The MathWorks, Inc.
parallel.internal.supportpackages.postInstallCmd(...
    parallel.internal.types.SupportPackageScheduler.Slurm);
end
