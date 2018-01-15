% SubmissionMode - Enumeration of submission modes for Generic profiles

%   Copyright 2017 The MathWorks, Inc.

classdef SubmissionMode < parallel.internal.types.NamedEnumeration
    enumeration
        Shared('shared', false, false)
        Remote('remote', true, false)
        Nonshared('nonshared', true, true)
    end
    
    properties (SetAccess = immutable)
        RequiresClusterHost
        RequiresRemoteJobStorageLocation
    end
    
    methods
        function obj = SubmissionMode(name, requiresClusterHost, requiresRemoteJobStorageLocation)
            obj@parallel.internal.types.NamedEnumeration(name);
            obj.RequiresClusterHost = requiresClusterHost;
            obj.RequiresRemoteJobStorageLocation = requiresRemoteJobStorageLocation;
        end
    end
    
    methods (Static)
        function enumVal = fromName(name)
            enumVal = parallel.internal.types.NamedEnumeration.getEnumFromName(...
                'parallel.internal.types.SubmissionMode', name);
        end
    end
end