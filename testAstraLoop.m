% Script which loops over different paramters and calls 'testAstra'


clear all;aclear;
%% Set parameter sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Experiment(1).name = 'Xeno4cell';
Experiment(1).subset = [12 7 2 5 11   6 1];
angleInc = [3 1 2 4 5 6];
horInc = [1];
filtLowFreq = [0];
maxIter = 250;
recType = {'SIRT_CUDA','CGLS_CUDA','SART_CUDA'};
maskDiskDecr = 8;
%% Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nHor = 1:numel(horInc)
    for nExp = 1:numel(Experiment)
        for nSub = 1:numel(Experiment(nExp).subset)
            for nAng = 1:numel(angleInc)
                for nLow = 1:numel(filtLowFreq)
                    for nRec = 1:numel(recType)
                        % Call function
                        testAstra('Experiment',Experiment(nExp).name,'SubsetIndex',Experiment(nExp).subset(nSub),'angleInc',angleInc(nAng),'horInc',horInc(nHor),'filtLowFreq',filtLowFreq(nLow),'maskDiskDecr',maskDiskDecr,'recType',recType{nRec},'maxIter',maxIter,'showFigs',0)
                    end
                end
            end
        end
    end
end