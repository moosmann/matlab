function c = nist_material_constants_elements( )
% Reads nist_material_constants_elements.txt and returns an cell array.
%
% RETURNS
% 1 x 6 cell array of cell arrays: { {atomic number Z}, {symbol}, {name},
% {Z/(mass A)}, {mean excitation energy I / (eV)}, {density / (g/cm3)}}

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [parpath, name] = fileparts( mfilename( 'fullpath' ) );
 filename = sprintf( '%s/%s.txt', parpath, name );
 fid = fopen( filename );
 c = textscan( fid, '%u%s%s%f%f%f', 'CommentStyle', '%', 'MultipleDelimsAsOne', 1);
 fclose( fid );
 