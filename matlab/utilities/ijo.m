function ijo(impath)

if nargin < 1
    system('imagej_opensequence');
else
    p = pwd;
    cd(impath)
    system('imagej_opensequence');
    cd(p)
end


end