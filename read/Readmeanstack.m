% ##		 READING ALL FILES FROM DIR TO THE SATCK
% ##
% ## 	MEANSTACK=readmeanstack(dir,nameprefix)
% ##
% ## 	where dir - directory where data is stored
% ## 	nameprefix - prefix of the file, 
% ##	f.e img_flat - automatically will read all files with this prefix
% ## 	all input parameters are STRING values
% ## 	Number of images in the folder displayed automatically
% ## 	output - 3D stack of images with nameprefix in the dir

function meanstack=readmeanstack(dir,nameprefix)

% Searching through the folder for all files which are matching with prefix
fname=sprintf('%s%s*',dir,nameprefix);
search=glob(fname);
fnames=strvcat(search);

[M,MM]=size(fnames); % number of files in directory
printf('Number of files in the folder=%d\n',M)

%load first image to get size
meanstack(:,:)=edfread(fnames([1],:));    

%load all other images to one file
for i=2:M
    meanstack(:,:)=meanstack(:,:)+edfread(fnames([i],:));
end

%normalize file
meanstack=meanstack(:,:)/M;

endfunction