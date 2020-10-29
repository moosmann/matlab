function RemoveTrailingSlash( string )
% Check for trailing '/' of input argument 'string' and append '/' if
% missing. The variable is changed implicitly in the caller function
% workspace using function assigin.

%% Check for trailing slash and append one when missing
if ~isempty( string )
    if string(end) == filesep
        string(end) = '';
    end
else
    return
end

assignin( 'caller', inputname(1), string)
