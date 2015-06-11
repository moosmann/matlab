%% Do <something> over a list of files obtained by 'ls' command.
%% It is exactly like pmedf_average4list(), just doing <something>
%% instead of averaging.
%%
%% That <something> here is: rotation by 90 degrees counter-clockwise,
%% i.e. rot90(a,3).
%%
%% Syntax:
%%	pmedf_FileList_rotate( outprefix, edfdir, edfnames )
%% where
%%	outprefix - can be a directory name (with slash)
%%	edfdir    - directory where are the input files
%%	edfnames  - input files, can include wildcards to be expanded by shell
%%
%% The list of files is obtained by the command
%%	unix( ['ls -1 ' edfdir '/' edfnames]
%% Name of the output file: the directory of the input files is stripped away
%% before applying the outprefix. Default of '' is current dir.
%%
%% Examples:
%%	pmedf_FileList_rotate( 'new/', '.', 'tomo_0???.edf' );
%%	pmedf_FileList_rotate( 'x',    '/data/id99/scan01/', '*.edf');
%%
%% Author: Petr Mikulik
%% Version: November 2003

function pmedf_FileList_rotate( outprefix, edfdir, edfnames )

if nargin ~= 3
    error('Usage: pmedf_FileList_rotate( outprefix, edfdir, edfnames )\n');
end

if strcmp(edfdir, '.')
    edfdir = [];
elseif length(edfdir) > 0
    if edfdir(length(edfdir))~='/' edfdir = [edfdir '/']; endif
end

if index(edfnames, ' ') > 0
    error('edfnames cannot contain space');
end

if length(edfdir) > 0
    command = ['ls -1 ', edfdir, edfnames];
else
    command = ['ls -1 ', edfnames];
end

if (0)
    % The original method (May 2002). Problem when there are too many files
    % and the command buffer gets overful -- no indication that not all files
    % were listed.
    [tmp, names] = unix(command);
else
    % Use a temporary file for redirecting 'ls' command to it, and then read
    % and delete that file.
    tmp = tmpnam();
    command = [command, '>', tmp];
    unix(command);
    fid = fopen(tmp, 'r');
    names = fscanf(fid, '%c');
    fclose(fid);
    unlink(tmp);
end

% proceed over all files
while (length(names) > 0)
    % set up current input name
    p = index(names, sprintf('\n'));
    if p==0
	name = names;
	names = [];
    else
	name = names(1:p-1);
        names = names(p+1:length(names));
    end
    p = rindex(name, '/');
    if p>0 name = name(p+1:length(name)); end
    inpname = name;
    if length(edfdir) > 0
	inpname = [edfdir inpname];
    end
    % set up current output file name
    outname = [outprefix name];

    if strcmp(inpname, outname)
	inpname
	outname
	error('FATAL ERROR IN INPUT: These files are the same!');
    end

    fprintf('Input:\t');
    [h, a] = pmedf_read(inpname);

    % Action: rotation
    % Get current header values
    d1 = pmedf_findInHeader(h, 'Dim_1', 'int');
    d2 = pmedf_findInHeader(h, 'Dim_2', 'int');
    p1 = pmedf_findInHeader(h, 'PSize_1', 'int');
    p2 = pmedf_findInHeader(h, 'PSize_2', 'int');
    % Write new header values
    newprefix = [outprefix, name];
    newprefix = newprefix(1:length(newprefix)-8);
    k = rindex(newprefix, '/');
    if (k>0) newprefix = newprefix(k+1:length(newprefix)); end
    newh = pmedf_putInHeader(h, 'PSize_1', p2);
    newh = pmedf_putInHeader(h, 'PSize_2', p1);
    newh = pmedf_putInHeader(h, 'Dim_1', d2);
    newh = pmedf_putInHeader(h, 'Dim_2', d1);
    newh = pmedf_putInHeader(h, 'prefix', newprefix);
    fprintf('  => New output: ');
    a = rot90(a,3);
    pmedf_write(outname, newh, a);
end

%eof pmedf_FileList_rotate.m
