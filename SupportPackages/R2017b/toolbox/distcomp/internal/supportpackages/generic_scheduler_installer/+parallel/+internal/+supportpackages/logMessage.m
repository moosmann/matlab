function logMessage(varargin)
%LOGMESSAGE - internal: logging for generic scheduler installer
% Call logMessage('showmessage', true) to enable logging.
% Call logMessage('showmessage', false) to disable logging.
% Call logMessage('something to log') to write the message 'something to log'.

% Copyright 2016 The MathWorks, Inc.

persistent showmessage;

% Initialize the persistent on first call.
% Logging is disabled by default.
if isempty(showmessage)
    showmessage = false;
end

if nargin == 2 && strcmp(varargin{1}, 'showmessage')
    validateattributes(varargin{2}, {'logical'}, {});
    showmessage = varargin{2};
    if showmessage
        iLogMessage('parallel.internal.supportpackage.logMessage() enabled');
    else
        iLogMessage('parallel.internal.supportpackage.logMessage() disabled');
    end
else
    if showmessage
        iLogMessage(varargin);
    end
end

end

function iLogMessage(varargin)
    disp(varargin{:});
end