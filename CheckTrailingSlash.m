function CheckTrailingSlash(String)
% Check for trailing '/' of input argument 'String' and append '/' if
% missing. The variable is changed implicitly in the caller function
% workspace using function assigin.

%% Check for trailing slash and append one when missing
if ~isempty(String)
    if String(end) ~= '/'
        String = [String '/'];
    end
else
    return
end

assignin('caller',inputname(1),String)
