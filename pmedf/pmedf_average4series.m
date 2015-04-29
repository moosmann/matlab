%% Average 2x2 over a whole series of files.
%%
%% Syntax:
%%	pmedf_average4series( outprefix, edfformat, numbers )
%% where
%%	outprefix - can be a directory name (with slash)
%%	edfformat - format string for input file names
%%	numbers	  - range of file name numbers
%%
%% If the input data contain a directory, then it is stripped away before
%% applying the outprefix.
%%
%% Examples:
%%	pmedf_average4series( 'new/', 'tomo%03i.edf', [0:800] );
%%	pmedf_average4series( 'x',    '/data/id99/scan01/x%03i.edf', [0:8]);
%%
%% Author: Petr Mikulik
%% Version: May 2002

function pmedf_average4series ( outprefix, edfformat, numbers )

if nargin ~= 3
  error('Usage: pmedf_average4series( outprefix, edfformat, numbers )\n');
  return
end

% proceed over all files
for k=1:length(numbers)
    % set up current input name
    name = sprintf(edfformat, numbers(k));
    nodirname = rindex(name, '/');
    nodirname = name(nodirname+1:length(name));
    % set up current output file name
    outname = [outprefix nodirname];

    fprintf('Input:\t');
    [h, a] = pmedf_read(name);
    h = pmedf_removeInHeader(h, 'prefix');

    % Action: average4
    [h4, a4] = pmedf_average4(h, a);
    [nr, nc] = size(a4);
    h4 = pmedf_putInHeader(h4, 'row_beg', '0', 16 );
    h4 = pmedf_putInHeader(h4, 'row_end', sprintf('%i',nr-1), 16 );
    h4 = pmedf_putInHeader(h4, 'col_beg', '0', 16 );
    h4 = pmedf_putInHeader(h4, 'col_end', sprintf('%i',nc-1), 16 );

    fprintf('  => Averaged Output: ');
    pmedf_write(outname, h4, a4);
end

%eof pmedf_average4series.m
