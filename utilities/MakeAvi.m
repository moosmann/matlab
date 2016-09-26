%aviobj = avifile('intensity.avi','compression','None');
%aviobj =
%avifile('intensity.m','compression','None','quality',100,'fps',2,'colormap',gray,'videoname','Intensity_pattern.m')

aviobj = avifile('int.msvc.avi','colormap',gray,'videoname','Intensity_pattern.m');
aviobj.FPS = 3;
aviobj.Quality = 100;
aviobj.COMPRESSION ='None';
for k=1:size(int,3)
    F      = int(:,:,k) ;
    aviobj = addframe(aviobj,F);
end;
aviobj = close(aviobj);