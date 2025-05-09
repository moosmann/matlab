function test_bang
    p = '/asap3/petra3/gpfs/p07/2024/data/11020415/scratch_cc/test';
    s = datestr(now, 'yyyymmddTHHMMSS');
    p2 = [p filesep s];
    CheckAndMakePath(p2)
    im = rand([1000,1000]);
    parfor n = 1:100
        fn = sprintf('%s/im_%04u.tif',p2,n);
        imwrite(im,fn,'tif');

    end
    cd(p2)

    %unix('ls');
    pwd
    unix('imagej_opensequence');

end