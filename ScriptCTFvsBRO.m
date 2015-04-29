
correction = 1;
blurring   = [8 8 1.5];
alpha      = 12;
pixelsize  = 1.1e-6;
%Ez         = {[30 .01] [10 .01] [10 .02] [10 .04] [10 .1] [10 1]};
Ez         = {[30 .01] [10 .01] [30 0.1] [10 .1] [10 .5] [10 1]};
%Ez         = {[10 .5]};

for phaseshift=[.01 .1 .5 1 6]
    %for phaseshift=[.01 .15  1 6]
    for x = [1:size(Ez,2)]
        energy   = Ez{x}(1);
        distance = Ez{x}(2);
        for counts = [0 20000 8000 2000 1000]
            %for counts = [0]
            [phase,int,phiCTF,phiBRO,ldp,namestring,intFT] = ...
            LenaCTF(counts,correction,blurring,alpha,phaseshift,[energy distance pixelsize]);
            assignin('base',['pha' namestring(1:5)],phase);
            assignin('base',['int',namestring],int);
            assignin('base',['intFT',namestring],intFT);
            assignin('base',['ctf',namestring],phiCTF);
            assignin('base',['bro',namestring],phiBRO);
        end;
    end;
end;
clear counts correction blurring alpha phaseshift ez energy distance pixelsize ...
    phase int phiCTF phiBRO ldp namestring x;
