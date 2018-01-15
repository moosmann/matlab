function registerMessageCatalog(installationFolder)
%REGISTERMESSAGECATALOG - internal: registers support package message catalog

% Copyright 2016 The MathWorks, Inc.

parallel.internal.supportpackages.logMessage(sprintf('registerMessageCatalog() resourceRoot %s', installationFolder));

% The message catalog has to be registered with MATLAB on a per-session basis.
% Re-registration of an existing catalog results in a harmless no-op.
matlab.internal.msgcat.setAdditionalResourceLocation(installationFolder);
end
