function byte = bytes_of_dtype( dtype )
% Returns the number of bytes of data type 'dtype'.
%
% byte = bytes_of_dtype( dtype )

% Convert data type to matlab data type
switch dtype
    case 'single' % float32, 'F'
        byte = 4;
    case 'double' % float64, 'D'
        byte = 8;
    case 'uint8'
        byte = 1;
    case 'uint16' %'U'
        byte = 2;
    case 'uint32'
        byte = 4;
    case 'uint64'
        byte = 8;
    otherwise
        error(['Data type not supported:' dtype]);
end