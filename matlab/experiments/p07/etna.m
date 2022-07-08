% beamtime-metadata-11012618
% "email": "Mattia.Pistone@uga.edu",
% "institute": "University of Georgia",
% "lastname": "Pistone",
% "proposalId": "20211140",
% "proposalType": "I",
% "title": "X-ray tomography of mineral-pore-glass contacts in volcanic rocks: unravelling the origin of paroxysms to improve their forecast at Mt etna"


par_path = '/asap3/petra3/gpfs/p07/2022/data/11012618/processed';
scan_subfolder = 'reco';
reco_subfolder = 'float_rawBin2';
stitched_volume_path = '';
% scan subsets as a cell array of strings
cc = {
{'etna_078_lp27h_d'
'etna_079_lp27h_e'}
{'etna_080_et21m_h'
'etna_081_et21m_i'}
{'etna_83_et28_forna_i'
'etna_084_et28_forna_j'}
{'etna_086_cse020421g_f'
'etna_087_cse020421g_g'}
{'etna_088_cse231021a_g'
'etna_089_cse231021a_h'}
{'etna_090_cse100222g_g'
'etna_091_cse100222g_h'}
{'etna_092_cse220222g_e'
'etna_093_cse220222g_f'}
{'etna_094_etn122ma_g'
'etna_095_etn122ma_h'}
{'ms_100_nex_3'
'ms_101_nex_3'}
};

for n = 2:numel(cc)
    c = cc{n};    
    fprintf( '\nheight scan %u', n)
    fprintf( '\n  %s', c{:} )
        
    c = cellfun(@(x) [par_path filesep x],c,'UniformOutput',false);
    stitched_volume_path = [c{end} '_stitched'];
    scan_path = c;
    stitch_volumes( scan_path, scan_subfolder, reco_subfolder, stitched_volume_path );
    
    
end
    
