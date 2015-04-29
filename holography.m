% Bilderserie mit Probe (Kunststoffkugeln mit 2 µm Durchmesser): http://dl.dropbox.com/u/1928462/Abstandsserie.zip
% Quelle bei gleichen Abständen: http://dl.dropbox.com/u/1928462/Quellen_Abstandsserie.zip
% 
% Die Abstände zwischen Quelle und Detektor sind in [mm] als erste 3 Ziffern des Dateinamens angegeben.
% Den Abstand zwischen Quelle und Probe muß ich am Montag nachschauen. Reiche diese Info also nach.
% 
% Die Aufnahmen wurden bei 260 eV gemacht. Als Quelle fungierte eine 400 nm große Lochblende.




clear all
SamDir = '/home/moosmann/holography/Abstandsserie';
SrcDir = '/home/moosmann/holography/Quellen_Abstandsserie';

dist = [366 369 372 375 376 396 416 426];
%sam = readstack(SamDir,'','tif',8,1,8,0.03);
%quelle = sam(:,:,end);%double(imread([SamDir '/quelle_100812_132151.tif']));
%src = readstack(SrcDir,'src','tif',8,1,8,0.03);
sam = readstack(SamDir,'','tif');
sam(:,:,end) = [];
src = readstack(SrcDir,'src','tif');
int = sam./src;

bf  = SineFilter([2048 2048],[0.26 0.366 0.7e-6],0.01);
out = Reco(int(:,:,1),2.5,[0.26 0.366 0.7e-6],[1 1]);
%itool(bf)
itool(out.ctfProjected)
itool(out.tieLO)
