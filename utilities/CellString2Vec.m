function vec = CellString2Vec( cc, imtype_str_flag )
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

format = cc{1}(end-2:end);

%% KIT camera: tiff
if strcmp( format, 'tif' )
    for nn = numel( cc ):-1:1
    if imtype_str_flag >= 0
            vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:5)) );                
        elseif imtype_str_flag == -1
                vec(nn) = str2double( cc{nn}(end-12:end-8) );                
        elseif imtype_str_flag == -2
                vec(nn) = str2double( cc{nn}(end-7:end-4) );                
        end        
    end
    
    %% KIT camera: raw
elseif strcmp( format, 'raw' )
    for nn = numel( cc ):-1:1
        if imtype_str_flag >= 0
            vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:5)) );                
        elseif imtype_str_flag == -1
                vec(nn) = str2double( cc{nn}(end-12:end-8) );                
        elseif imtype_str_flag == -2
                vec(nn) = str2double( cc{nn}(end-7:end-4) );                
        end        
    end
    %% EHD camera
else
    for nn = numel( cc ):-1:1
        vec(nn) = str2double( cc{nn}(end-8:end-4) );
        
    end
end
