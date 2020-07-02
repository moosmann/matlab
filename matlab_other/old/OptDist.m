%Script for quick phase retrieval
clear all
%% Parameters for reading the data.
DataPath = '/mntdirect/_data_visitor/mi1057/bm05/pc/opt_xeno_stage12_F_/';
ScanPrefix = 'opt_xeno_stage12';
HPFT = 0.03;

NumExp = numel(dir([DataPath ScanPrefix '*dist00_exp*' ]));
NumDist = numel(dir([DataPath ScanPrefix '*dist*_exp00*' ]));
fprintf('\nNumber of exposure times: %u\n',NumExp)
fprintf('Number of distances: %u\n',NumDist)

for ii = NumDist:-1:1
    arDat = dir(sprintf([DataPath ScanPrefix '*dist%02u_exp*' ],ii-1));
    for jj = NumExp:-1:1
         % Read darks.
        arDark = dir(sprintf([DataPath 'dark_' ScanPrefix '*_exp%02u*' ],jj-1));
        if ii == NumDist
            for rr = numel(arDark):-1:1
                dark(:,:,rr) = pmedfread(arDark(rr).name);
            end
            darkMedian(:,:,jj) = FilterHotPixel(median(dark,3),HPFT);
        end
        % Read refs.
        arRef = dir(sprintf([DataPath 'ref_' ScanPrefix '*dist%02u_exp%02u*' ],ii-1,jj-1));
        for rr = numel(arRef):-1:1
            fprintf('%02u  %02u  %02u   ',ii,jj,rr)
            ref(:,:,rr) = pmedfread(arRef(rr).name);
        end
        fprintf('\n')
        refMedian = FilterHotPixel(median(ref,3)-darkMedian(:,:,jj),HPFT);
        % Read data.
        data = pmedfread(sprintf([DataPath arDat(jj).name]));
        im{ii,jj} = FilterHotPixel((data-darkMedian(:,:,jj))./refMedian,HPFT);
        edfwrite(sprintf([DataPath 'ffcor/int_dist%02u_exp%02u.edf'],ii,jj),im{ii,jj},'float32');
    end
end

