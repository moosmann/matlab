for ii =[21 18 15 12 9];
    edfwrite(sprintf('XenopusStage25Agar_phaseCTFprojected_alpha%03.2f.edf',out(ii).alpha),out(ii).ctfProjected,'float32');
end