function vec = CellString2Vec( cc  )
% Convert cell of string which contain numeric expressions into an array of
% the corresponding numeric values (double)
%
% INPUT
% cc: cell. cell of string
%
% OUTPUT
% vec: array of double
%
% Written by Julian Moosmann
% vec = CellString2Vec(cc)
len = [];
format = cc{1}(end-2:end);


% position of running index
re = regexp( cc{1}, '\d{6,6}');
len = 6;
if isempty( re )
    re = regexp( cc{1}, '\d{5,5}');
    len = 5;
end
if isempty( re )
    re = regexp( cc{1}, '\d{4,4}');
    len = 4;
end
if isempty( re)
    len = [];
end
if numel( re ) >= 1
    re = re(end);
    imtype_str_flag = re;
elseif strcmpi( cc{1}(end-6:end-4), 'ref' )
    imtype_str_flag = '64'; % 0
elseif strcmpi( cc{1}(end-11:end-9), 'ref' )
    imtype_str_flag = '119'; % 1
end

%% KIT camera: tiff
if strcmp( format, 'tif' )
    for nn = numel( cc ):-1:1
        if imtype_str_flag >= 0
            if ~isempty( len )
                vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:len-1)) );
            else
                vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:5)) );
            end
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
            if ~isempty( len )
                vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:len-1)) );
            else
                vec(nn) = str2double( cc{nn}(imtype_str_flag + (0:5)) );
            end
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
