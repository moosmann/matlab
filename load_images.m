% Function to load medipix (txt) images to matlab
%
%  function [n, A] = get_mpx_pictures ( directory, fpattern , filetype )
%
% Output:
%  n - number of loaded images
%  A - images, A(m,n,i)  - i is index of image m,n - number of pixels
%                                           horizontaly and verticaly
%  Input:
% directory - path to the files including last slash
% filepattern  - file pattern 
% filetype: 'mpx'   - medipix images
%           'tif'   - tif images  
%           'edf'   - edf images
%
% Written by P. Vagovic, 23. 7. 2010


function [n, A] = load_images ( directory, filepattern , filetype )

  % Check if last character of 'folder' string is a slash (/).
  if (directory(end)~='/'), directory = [directory,'/'];end

  datafiles= dir([directory filepattern]);
  n = length(datafiles);
  
  %medipix images, ASCCI matrix of in tensity
  if(filetype=='mpx') 
  
    for i=1:n
      A(:,:,i) = importdata([directory datafiles(i).name]);
      disp(['Reading image: ' directory datafiles(i).name ]);
    end
 
  %standard tif file from for example PCO4000
  elseif(filetype=='tif') 
     
     A = imread([directory datafiles(1).name]);
     [dim1,dim2] = size(A);
     A = zeros(dim1,dim2,n);  

    for i=1:n
      A(:,:,i) = imread([directory datafiles(i).name]);
      disp(['Reading image: ' directory datafiles(i).name ]);
    end
  
   %EDF 
   elseif(filetype=='edf') 
        for i=1:n
            [h, A(:,:,i)] = pmedf_read([directory datafiles(i).name]);
      
        end
   else
        disp(['Error: Filetype ' filetype ' is not supported!!!']); 
   end
   
   return
