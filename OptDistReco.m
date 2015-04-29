clear all
%% Parameters
alpha = [5 2.5];
evalTIElo = 1;
evalTIEpnlo = 0;
evalCTF = 0;
BinaryFilterThreshold = 0.01;
offset = 0.023;
xc_farest = 0.400;
xc_closest = 0.050;
dxc = (xc_farest-xc_closest)/7;

DataPath = '/mntdirect/_data_visitor/mi1057/bm05/pc/opt_xeno_stage12_F_/';
fprintf('Processing:')
for ii = 8:-1:1
    for jj = 7:-1:1
        fprintf('%u %u   ',ii,jj)
        im0 = pmedfread(sprintf('%sffcor/int_dist%02u_exp%02u.edf',DataPath,ii,jj));
        im = im0(300:1748,300:1748);
        out = Reco(im,alpha,[22 (offset+xc_farest-ii*dxc) 0.75e-6],evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold);
        pctf = out.ctfProjected;
        tie = out.tieLO;
        edfwrite(sprintf([DataPath 'phase/pctf_dist%02u_exp%02u.edf'],ii,jj),pctf,'float32');
        edfwrite(sprintf([DataPath 'phase/tieLO_dist%02u_exp%02u.edf'],ii,jj),tie,'float32');
        %edfwrite(sprintf([DataPath 'phase/pctf_LFR_dist%02u_exp%02u.edf'],ii,jj),RemoveLowFreq(pctf),'float32');
    end
    fprintf('\n')
end
