%% edfs2vtk.m
%% Convert a series of edf files into a volume vtk file.
%% The input files are given as a list of strings, i.e. files{1}, files{2},
%% etc., see examples below.
%%
%% Syntax:
%%	edfs2vtk(out_vtk_filename, edf_filelist [, opts]);
%% where the opts structure can have the following fields (they are processed
%% in the given order):
%%	opts.roi = [{0,0 | xfrom, xto}; {0,0 | yfrom, yto}]
%%				    ... do (not)cut subimage
%%					negative indices are wrt the opposite
%%					corner, which is -1,-1
%%	opts.binning = [binx, biny] ... do averaging along x,y pixel regions
%%	opts.bin_z = binz	    ... do averaging over subsequent images
%%	opts.zrange = [zfrom,zto]   ... set all voxel values into this range
%%	opts.datatype = 'type'	    ... output datatype: uchar, ushort, float
%%	opts.zspan = {1|[zfrom,zto]}... let zrange span this range of datatype
%%
%% Notes:
%%   1. The binning and bin_z options are really averaging, not binning.
%%   2. Using opts.zrange makes sense for setting low background (or intense
%%	peaks) to a unique level of zfrom (or zto) values.
%%   3. You can choose output datatype different than that in edf files;
%%	supported are uchar, ushort, float.
%%   4. If the output datatype has shorter range of values than that in edf
%%	files (e.g. conversion ushort => uchar), then you can span (rescale)
%%	the zrange=[zfrom,zto] to output range zspan:
%%	    - opts.zspan=0 or it is undefined: no rescaling
%%	    - opts.zspan=1: span over the whole range of opts.datatype, e.g.
%%			    [0,255] for uchar
%%	    - opts.zspan=[0,5]: rescales the zrange into range [0,5].
%%
%% Examples:
%%	edfs2vtk('out1.vtk', wildcard2filelist('volume/slice???.edf'));
%%	edfs2vtk('out2.vtk', wildcard2filelist('data/scan_0??[0,5].edf'));
%%
%%	files = wildcard2filelist('data/sample_0deg_*.edf');
%%	edfs2vtk('out3.vtk', files(1:4:end));
%%	clear opts;
%%	opts.binning=[6,6]; opts.bin_z=2;
%%	opts.roi=[256,-256;0,0]; % for 1024x1024, it means the central line
%%	edfs2vtk('out4.vtk', files, opts);
%%
%%	opts.datatype='uchar'; opts.zrange=[100,300];
%%	opts.zspan=1;  or  opts.zspan=[0,200];
%%
%% Author: Petr Mikulik
%% Version: February 2005
%% History:
%%	18. 2. 2005: Added opts.datatype, opts.zrange, opts.zspan
%%	15. 2. 2005: Added opts.roi and opts.bin_z.
%%      23. 9. 2004: Original version.

function edfs2vtk ( out_vtk_filename, filelist, opts )

if nargin<2 || nargin>3
    error('Usage: edfs2vtk(out_vtk_filename, edf_filelist [, opts]');
end

if ~isstr(out_vtk_filename)
    error('Output filename must be a string');
end

if ~iscell(filelist) | length(filelist)==0
    error('Wrong filelist (string cell array required)');
end

% Default values
binning=[1,1];
bin_z=1;
roi_x=[];
roi_y=[];
outdatatype_m=[];
zrange=[];
zspan=[];

% Parse options
if nargin==3
    if ~isstruct(opts)
	error('Options must be a structure');
    end
    a=fieldnames(opts);
    for k=1:length(a)
	if strcmp(a{k},'binning')
	    binning=opts.binning;
	    if ~isvector(binning) | length(binning)~=2
		error('opts.binning must be vector with two positive integers');
	    end
	elseif strcmp(a{k},'bin_z')
	    bin_z=opts.bin_z;
	    if ~isscalar(bin_z) | bin_z<0
		error('opts.bin_z must be a positive integer');
	    end
	    if bin_z==0 bin_z=1; end % accept also 0 as no binning
	elseif strcmp(a{k},'roi')
	    if any(size(opts.roi)~=[2,2])
		error('opts.roi must be a 2x2 matrix');
	    end
	    if all(opts.roi(1,:)~=0) roi_x=opts.roi(1,:); end
	    if all(opts.roi(2,:)~=0) roi_y=opts.roi(2,:); end
	elseif strcmp(a{k},'datatype')
	    outdatatype_m=pmvtk_getDatatype(opts.datatype);
	    outdatatype_m=outdatatype_m.m;
	elseif strcmp(a{k},'zrange')
	    zrange=opts.zrange;
	    if ~isvector(zrange) | length(zrange)~=2
		error('opts.zrange must be a vector with just two elements');
	    end
	elseif strcmp(a{k},'zspan')
	    zspan=opts.zspan;
	    if length(zspan)==1
		if (zspan~=0) zspan=1;
		else zspan=[];
		end
	    elseif length(zspan)~=2
		error('opts.zspan must be 0, 1, or [from,to]');
	    end
	else
	    error(['Unknown option "opt.', a{k}, '"']);
	end
    end
end
if ~isempty(zspan) & isempty(zrange)
    error('If you define opts.zspan, you must define opts.zrange as well');
end

% Read the 1st file
[h, a] = pmedf_read(filelist{1});
% Update roi
[dim1, dim2] = size(a);
if ~isempty(roi_x)
    if (roi_x(1) < 0) roi_x(1)=roi_x(1)+dim1+1; end % negative roi = wrt opposite corner
    if (roi_x(2) < 0) roi_x(2)=roi_x(2)+dim1+1; end
    if (roi_x(2) > dim1) roi_x(2)=dim1; end % protection against too big range
    if (roi_x(1) > roi_x(2)) roi_x(1)=1; end
end
if ~isempty(roi_y)
    if (roi_y(1) < 0) roi_y(1)=roi_y(1)+dim2+1; end % negative roi = wrt opposite corner
    if (roi_y(2) < 0) roi_y(2)=roi_y(2)+dim2+1; end
    if (roi_y(2) > dim2) roi_y(2)=dim2; end % protection against too big range
    if (roi_y(1) > roi_y(2)) roi_y(1)=1; end
end
% Cut the subimage
if ~isempty(roi_x) & ~isempty(roi_y) a=a(roi_x(1):roi_x(2),roi_y(1):roi_y(2));
elseif ~isempty(roi_x) a=a(roi_x(1):roi_x(2),:);
elseif ~isempty(roi_y) a=a(:,roi_y(1):roi_y(2));
end

% Make x- and y-axes.
[dim1, dim2] = size(a);
[newdim1, newdim2] = size(a);
psize1 = pmedf_findInHeader(h, 'PSize_1', 'float');
psize2 = pmedf_findInHeader(h, 'PSize_2', 'float');
if isempty(psize1) psize1=1; end
if isempty(psize2) psize2=1; end
if binning(1) ~= 1
    psize1=psize1*binning(1);
    newdim1=floor(dim1/binning(1));
    dim1=binning(1)*newdim1;
end
if binning(2) ~= 1
    psize2=psize2*binning(2);
    newdim2=floor(dim2/binning(2));
    dim2=binning(2)*newdim2;
end
%x_axis = psize1*(0:newdim1-1);
%y_axis = psize2*(0:newdim2-1);

if (any(binning>1))
    a = matbinning(a, binning) / prod(binning);
end

% Number of input files and output slices
nfiles = length(filelist);
% slices to write = nfiles over bin_z chunks
if bin_z>1
    nslices=floor(nfiles/bin_z);
    nfiles=bin_z*nslices;
else
    nslices=nfiles;
end

% Make z-axis.
% We should extract the rtn motor instead by pmedf_findInHeader().
% However, by default MayaVi uses true aspect ratio of volume axes, I don't
% know how to change the scaling to "fit to window", so we should keep it
% in the same ratio as spatial axes.
% If someone knows how to scale axes within MayaVi, let me know!
% z = 0:nslices-1;
%z_axis = linspace(0, 0.5*(x(end)+y(end)), nslices);

indatatype = pmedf_findInHeader(h, 'DataType', 'string');

% remember input datatype
switch indatatype
    case {'FloatValue', 'Float', 'Float32'}, indatatype_m='float32'; indatatype_vtk='float';
    case 'UnsignedShort', indatatype_m='ushort'; indatatype_vtk='unsigned_short';
    case 'UnsignedChar', indatatype_m='uchar'; indatatype_vtk='unsigned_char';
    otherwise error(['unimplemented datatype "',indatatype,'", sorry']);
end

% set output datatype
if isempty(outdatatype_m) outdatatype_m=indatatype_m; end
switch outdatatype_m
    case 'float32', outdatatype_vtk='float';
		    if length(zspan)==1 zspan=[]; end
    case 'ushort', outdatatype_vtk='unsigned_short';
		    if length(zspan)==1 zspan=[0,65535]; end
    case 'uchar', outdatatype_vtk='unsigned_char';
		    if length(zspan)==1 zspan=[0,255]; end
    otherwise error(['unimplemented output datatype "',outdatatype_m,'", sorry']);
end

% Open the output file
fprintf('==> Writing VTK file "%s"', out_vtk_filename);
if (any([binning,bin_z]>1))
    fprintf('\n==> with (%ix%i)x%i binning (thus "%s" volume is %ix%ix%i)', binning,bin_z, outdatatype_m, newdim1,newdim2,nslices);
end
fprintf('\n');

f = fopen(out_vtk_filename, 'wb');

header = [...
'# vtk DataFile Version 2.0\n',...
'GID volume\n',...
'BINARY\n',...
'DATASET STRUCTURED_POINTS\n',...
'DIMENSIONS %i %i %i\n',...
'SPACING %g %g %g\n',...
'ORIGIN 0 0 0\n',...
'POINT_DATA %lu\n',...
'SCALARS scalars %s 1\n',...
'LOOKUP_TABLE default\n'...
];

k=newdim1*newdim2*nslices;
fprintf(f, header, newdim1,newdim2,nslices, psize1,psize2,1, k, outdatatype_vtk);

% Read edf files, write slices
mymin=zeros(1,nslices); % allocation
mymax=zeros(1,nslices);
s=0; % slice number (1..nslices)
ks=0; % file number in the slice (1..bin_z)
for k=1:nfiles
    if k>1 % first file of the series has been already read
	[h,a]=pmedf_read(filelist{k});
	% Cut the subimage
	if ~isempty(roi_x) & ~isempty(roi_y) a=a(roi_x(1):roi_x(2),roi_y(1):roi_y(2));
	elseif ~isempty(roi_x) a=a(roi_x(1):roi_x(2),:);
	elseif ~isempty(roi_y) a=a(:,roi_y(1):roi_y(2));
	end
	if any(binning>1)
	    a = matbinning(a, binning) / prod(binning);
	end
    end

    if bin_z>1 % continue reading if binning over z
	ks=ks+1;
	if ks==1
	    asum=a;
	else
	    asum=asum+a; % Matlab does not have +=
	end
	if ks<bin_z % loop read next file
	    continue;
	end
	ks=0; % bin_z files read
	a=asum/bin_z;
    end

    % change voxel values into zrange
    if ~isempty(zrange)
	a(a<zrange(1))=zrange(1); % set small values to zrange(1)
	a(a>zrange(2))=zrange(2); % set large values to zrange(2)
    end

    % span zrange into the full uchar/ushort range
    if ~isempty(zspan)
	a=zspan(1)+(a-zrange(1))/(zrange(2)-zrange(1))*(zspan(2)-zspan(1));
    end

    s=s+1;
    count1 = fwrite(f, a, outdatatype_m, 0, 'ieee-be'); % VTK requires big endian
    mymin(s) = min(a(:));
    mymax(s) = max(a(:));
    fprintf('==> Slice %2i: min voxel value is %g, max is %g\n', s, mymin(s), mymax(s));
end

fclose(f);

fprintf('\nConversion  "%s"\n\t... "%s"\n\t==> "%s"\nhas been done.', filelist{1}, filelist{nfiles}, out_vtk_filename);
fprintf('\nVolume is %ix%ix%i of "%s"', newdim1, newdim2, nslices, outdatatype_m);
fprintf(', min voxel value is %g, max is %g.\n', min(mymin), max(mymax));

% eof edfs2vtk.m

