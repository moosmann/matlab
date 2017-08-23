ca
clear p v d im
p = single(Readstack('.',1,'proj_*_1.tif'));


for nn = size(p,3):-1:1
    v(:,:,nn) = Binning( p(:,:,nn), 2 );
    
end


v = normat(v(51:400,501:1050,:));
%nimplay(v)
save_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/videos/synchroloadMeeting/11001978/radio_03/';
for nn = 1:size(v,3)    
    filename = sprintf( '%sradio_%06u.tif', save_path, nn);
    imwrite( uint8( (2^8 - 1) * v( :, :, nn) ), filename );

end

for nn = size(p,3):-1:2
    im = abs(v(:,:,nn) - v(:,:,nn-1));
    domain(im)
    d(:,:,nn -1) = im;
end
%nimplay(v)