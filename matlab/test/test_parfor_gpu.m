OpenParpool(20);
% if ~exist('vol','var')
%     vol = ones([5000,5000,200]);
% end
im0 = rand([5000 5000]);
for numWorkers = 1:20
    fprintf('\n# Workers: %u\n',numWorkers)
    tic;
    parfor (nn = 1:40,numWorkers)
        im = im0;
        a = gpuArray(im);
        a = FilterPixel(a);
        im = gather(a);
    end
    toc
end

fprintf('\nEND')