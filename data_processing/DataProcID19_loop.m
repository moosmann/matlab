
PhaseMethod = {'qp'};
%PhaseMethod = {'qp','qp2','ctf','tie'};
RegPar = [2 2.5];
filterSinoAfterPhase = [1 2 3];

for pm=1:numel(PhaseMethod)
    for rp = 1:numel(RegPar)
        for fsap = 1:numel(filterSinoAfterPhase)
            DataProcID19('PhaseMethod',PhaseMethod{pm},'RegPar',RegPar(rp),'filterSinoAfterPhaseRetrieval','filterSinoAfterPhase',filterSinoAfterPhase(fsap))
            %celldisp({'PhaseMethod',PhaseMethod{pm},'RegPar',RegPar(rp)})
        end
    end
end