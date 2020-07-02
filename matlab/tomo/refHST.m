function ref = refHST(folder,ref_list,write_edf)
% Compute refHST images. Takes median of the 15 ref images taken at
% 0000,0500,1000 adn 1500.

% folder - string: path to ref files.
% ref_list - list of integers: image number where the ref images were taken.
% write_edf - boolean: write refHST*.edf to folder.
    
    if nargin<1,folder='.';end;
    if nargin<2,ref_list = [0:500:1500];end;
    if nargin<3,write_edf=1,end;
    
% Check if last character of  string 'folder' is a slash (/).
if (folder(end)~='/'), folder = [folder,'/'];end

dimx = 2048;
dimy = 2048;
number_of_ref_files = 15;

% Loop over ref files.
for jj=ref_list,
     ref=zeros(dimx,dimy,number_of_ref_files);
     for ii=0:14,
         ref(:,:,ii+1)=edfread([folder 'ref' num2str(ii,'%04i') '_' num2str(jj,'%04i') ...
                             '.edf']);
     end;
     ref = median(ref,3);
     if write_edf,
         edfwrite([folder 'refHST' num2str(jj,'%04i') '.edf'],ref, ...
                  'float32');
     end;
end;
