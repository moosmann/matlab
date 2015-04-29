function m = Delta_Dot(res,dot_size,dot_magnitude)
% Creates grid of zeros with dimensions res=[dim_x,dim_y] and a central dot
% in the middle of size and magnitude dot_size and dot_magnitude.
    
    if (nargin<1) || isempty(res), res = [512,512]; end;
    if (nargin<2) || isempty(dot_size), dot_size = 1; end;
    if (nargin<3) || isempty(dot_magnitude), dot_magnitude = res(1); end;
    
    if (length(res)==1), res = [res,res]; end;  
    m = zeros(res);
    nx = res(1);
    ny = res(2);
    [x,y] = meshgrid(-floor(ny/2):floor(ny/2)-1,-floor(nx/2):floor(nx/2)-1);
    % [X,Y] = MESHGRID(x,y): The rows of the output array X are copies of the
    % vector x and the columns of the output array Y are copies of the
    % vector y.
    r = sqrt(x.^2 + y.^2);
    m(find(r<dot_size)) = dot_magnitude;
    
    
    
        
    