function SimulationSiemensStarVol(z,res,write_phex)

    if nargin<3,write_phex=0;end;
    
e=30;
%z=0.3;
%res=512;
d=1e-7;
s=1.5;
h=10;

prefix = sprintf(['/home/moosmann/data/SiemensStar/rawfiles/' ...
                  'E%gkev_z%03gcm_res%04g_delta%04gnm_10sigma%g_hsize%02g'],e,100*z,res,10^9*d,10*s,h);
prefix2 = sprintf(['/home/moosmann/data/SiemensStar/rawfiles/' ...
                  'E%gkev_z%03gcm_res%04g_delta%04gnm_10sigma%g_hsize%02g'],e,100*0,res,10^9*d,10*s,h);

fprintf('%s\n',prefix);
[phex,int,phi,partab] = SimulationSiemensStar(e,z,res,1,h,s,32,d);
[dimx dimy] = size(phex);

mexVolWrite(sprintf('%s_phase_reco.vol',prefix),phi,'float32');
mexVolWrite(sprintf('%s_intensity.vol',prefix),reshape(int,[1 dimx dimy]),'float32');
if write_phex==1,
    mexVolWrite(sprintf('%sphase_exact.vol',prefix2),reshape(phex,[1 dimx ...
                        dimy]),'float32');
end;
fid = fopen(sprintf('%s_parameter_file.txt',prefix),'wt');
fprintf(fid,sprintf(['                                 resolution    ' ...
                    'pixelsize                              exact phase             ' ...
                    '' ']reconstructed phase      real space error\n  E    ' ...
                    'z   detsiz    delta spok  binnd unbind  binned unbinnd ' ...
                    'alph sigm hsiz binf     min     max   range     min     ' ...
                    'max   range       <bro>    <brocor>\n%s'],partab));
fclose(fid);