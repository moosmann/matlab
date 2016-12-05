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
% Written by Julian Moosmann, 2016-11-08. Modified:
%
% vec = CellString2Vec(cc)

%cellpos = regexp( cc, '\d' );

for nn = numel( cc ):-1:1
   %vec(nn) = str2double( cc{nn}(cellpos{nn}) );
   vec(nn) = str2double( cc{nn}(end-8:end-4) );
   
end