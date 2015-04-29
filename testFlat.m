% flatStr = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/data/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/proj_01400.tif';
% flat = double(imread(flatStr));
% flat = FilterHotPixel(flat,0.03);
[dim1 dim2]=size(flat);
f=flat-mean(flat(:));
fp=padarray(f,0,'symmetric');
[dimx dimy]=size(f);
figure(1), imshow2(fp');
fpf = fft(fp,[],1);
%fpf = fftshift(fpf);
domain(abs(fpf),'abs(fpf)');
figure(2), imshow2(abs(fpf)');
%fpf2 = fpf;
sz = 70;
%fpf2([1:sz (2*dimx-sz):2*dimx],:)=0;
fpf2 = fpf(sz:(2*dimx-sz),:);
figure(3), imshow2(abs(fpf2)')
f2=ifft(fpf2,[],1);
domain(real(f2))
domain(imag(f2))
figure(4), imshow2(real(f2'))