n = 2000;

fprintf('\n host: %s',getenv('HOSTNAME'))
fprintf('\n volume size: %u^3',n)
if ~exist('proj','var') || size(proj,1)~=n
    fprintf('\n create volume: ')
    t = toc;
    proj = 10*ones([n n n],'single');
    proj_shape1 = size(proj,1);
    fprintf('%9.6f s',toc - t)

    padding_method = 'symmetric';
    filt = ones([2*proj_shape1, 1],'single');
end
t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;
t5 = 0;
t6 = 0;
t7 = 0;
t8 = 0;

fprintf('\n\nSERIAL: ')
for nn = 1:min([1000,size(proj,1)])
    % if mod(nn,10)==0
    %     fprintf(' %u',nn)
    % end
    % fprintf(' %u',nn)
    % if mod(nn,60)==0
    %     fprintf(' \n')
    % end
    t = toc;
    im = proj(:,:,nn);
    t1 = t1 + toc - t;

    t = toc;
    im = padarray(im,1 * [proj_shape1 0 0],padding_method,'post');
    t2 = t2 + toc - t;

    t = toc;
    im = fft(im,[],1);
    t3 = t3 + toc - t;

    t = toc;
    im = fun_times(im,filt);
    t4 = t4 + toc - t;

    t = toc;
    im = ifft(im,[],1,'symmetric');
    t5 = t5 + toc - t;

    t = toc;
    im = real(im);
    t6 = t6 + toc - t;

    t = toc;
    im = im(1:proj_shape1,:,:);
    t7 = t7 + toc - t;

    t = toc;
    proj(:,:,nn) = im;
    t8 = t8 + toc - t;
end

fprintf('\n %20s: %9.6f','assign image',t1)
fprintf('\n %20s: %9.6f','padding',t2)
fprintf('\n %20s: %9.6f','fft',t3)
fprintf('\n %20s: %9.6f','fun_times',t4)
fprintf('\n %20s: %9.6f','ifft',t5)
fprintf('\n %20s: %9.6f','real',t6)
fprintf('\n %20s: %9.6f','crop',t7)
fprintf('\n %20s: %9.6f','reassign image',t8)
fprintf('\n %20s: %9.6f','total',t1+t2+t3+t4+t5+t6+t7+t8)

%% Parallel
fprintf('\n\nPARALLEL')

OpenParpool(50);
t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;
t5 = 0;
t6 = 0;
t7 = 0;
t8 = 0;tic;
parfor nn = 1:min([100,size(proj,1)])

    t = toc;
    im = proj(:,:,nn);
    t1 = -t + toc + t1;

    t = toc;
    im = padarray(im,1 * [proj_shape1 0 0],padding_method,'post');
    t2 = -t + toc + t2;

    t = toc;
    im = fft(im,[],1);
    t3 = -t + toc + t3;

    t = toc;
    im = fun_times(im,filt);
    t4 = -t + toc + t4;

    t = toc;
    im = ifft(im,[],1,'symmetric');
    t5 = -t + toc + t5;

    t = toc;
    im = real(im);
    t6 = -t + toc + t6;

    t = toc;
    im = im(1:proj_shape1,:,:);
    t7 = -t + toc + t7;

    t = toc;
    proj(:,:,nn) = im;
    t8 = -t + toc + t8;
end
t = toc;
fprintf('\n %20s: %9.6f','assign image',t1)
fprintf('\n %20s: %9.6f','padding',t2)
fprintf('\n %20s: %9.6f','fft',t3)
fprintf('\n %20s: %9.6f','fun_times',t4)
fprintf('\n %20s: %9.6f','ifft',t5)
fprintf('\n %20s: %9.6f','real',t6)
fprintf('\n %20s: %9.6f','crop',t7)
fprintf('\n %20s: %9.6f','reassign image',t8)
fprintf('\n %20s: %9.6f','total sum',t1+t2+t3+t4+t5+t6+t7+t8)
fprintf('\n %20s: %9.6f','total',t)

% %% Parallel 
% fprintf('\n\nPARALLEL loop with increasing numbers of workers')
% num_gpu_max = 50;
% tic;
% tngpu = zeros([num_gpu_max,1]);
% for ngpu = 1:num_gpu_max
%     t = toc;
%     parfor (nn = 1:min([100,size(proj,1)]),ngpu)
%         im = proj(:,:,nn);
%         im = padarray(im,1 * [proj_shape1 0 0],padding_method,'post');
%         im = fft(im,[],1);
%         im = fun_times(im,filt);
%         im = ifft(im,[],1,'symmetric');
%         im = real(im);
%         im = im(1:proj_shape1,:,:);
%         proj(:,:,nn) = im;
%     end
%     tngpu(ngpu) = toc - t;
%     fprintf('\n %3u: %9.6f',ngpu,tngpu(ngpu))
% end
% fprintf('\n')
% 
% figure('Name','computation time vs poolsize')
% plot(1:num_gpu_max,tngpu)
% xlabel('poolsize')
% ylabel('compuation time / s')

fprintf('\n')

    
