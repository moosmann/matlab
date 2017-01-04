function vec = CellString2Vec(cc)
% Convert cell of string which contain numeric expressions into an array of
% the corresponding numeric values (double)
% 
% INPUT
% cc: cell. cell of string
%
% OUTPUT
% vec: array of double
%
% Written by Julian Moosmann, 2017-01-03. Modified:
%
% vec = CellString2Vec(cc)

%% KIT camera
if strcmp( cc{1}(end-2:end), 'tif' )
    for nn = numel( cc ):-1:1
    vec(nn) = str2double( cc{nn}(end-7:end-4) );   
    end
    
%% EHD camera    
else    
    for nn = numel( cc ):-1:1   
    vec(nn) = str2double( cc{nn}(end-8:end-4) );
   
    end
end
