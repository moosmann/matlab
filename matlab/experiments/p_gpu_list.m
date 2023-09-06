function gpu_index_list = p_gpu_list(sino, tomo, par, verbose)
% Return list of indices of GPU devices to be used in a parpool loop. The
% list takes care of the available memory per GPU device.

if nargin < 4
    verbose = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mem_avail_gpu = par.mem_avail_gpu;
mem_total_gpu = par.mem_total_gpu;
gpu_index = par.gpu_index;
poolsize_gpu_limit_factor = par.poolsize_gpu_limit_factor;
vol_shape = tomo.vol_shape;

% GPUs used
if isempty(gpu_index)
    gpu_index = 1:gpuDeviceCount;
end
% CPU pool size
num_gpu_used = numel(gpu_index);
p = gcp('nocreate');

parpool_size = p.NumWorkers;

% Acquire GPU memory consumption
if isempty(mem_avail_gpu) || numel(mem_avail_gpu) ~= num_gpu_used
    fprintf('\n Acquiring GPU memory consumption')
    for mm = num_gpu_used:-1:1
        nn = gpu_index(mm);
        gpu = gpuDevice(nn);
        gpu.reset;
        mem_avail_gpu(nn) = gpu.AvailableMemory;
        mem_total_gpu(nn) = gpu.TotalMemory;        
    end    
end

% Create list of GPU indices using all GPUs and accounting for available
% memory
gpu_mem_requ_per_reco = 1.05 * (vol_shape(1)*vol_shape(2) + numel(sino))*4;
recos_per_gpu = floor( poolsize_gpu_limit_factor  * mem_avail_gpu / gpu_mem_requ_per_reco);
gpu_index_list = zeros([1, parpool_size]);
tmp = recos_per_gpu;
for n = 1:parpool_size    
%     [v,p] = max(tmp);
%     if v == 0
%         break;
%     end
%     tmp(p) = tmp(p) - 1;
%     gpu_index_list(n) = gpu_index(p);
    m = mod(n - 1,num_gpu_used) + 1;
    if tmp(m) > 0
        gpu_index_list(n) = gpu_index(m);
        tmp(m) = tmp(m) - 1;
    end
end
gpu_index_list(gpu_index_list==0) = [];

if verbose
    fprintf('\n Parpool size: %u', parpool_size)
    fprintf( '\n GPUs to use:')
    fprintf( ' %u', gpu_index)
     for nn = 1:num_gpu_used
        ma = mem_avail_gpu(nn)/1024^3;
        mt = mem_total_gpu(nn)/1024^3;
        r = 100 * ma / mt;
        fprintf( '\n GPU %u: memory: total: %.3g GiB, available: %.3g GiB (%.2f%%)', gpu_index(nn), mt, ma, r)
    end
    fprintf('\n GPU poolsize limit factor : %g', poolsize_gpu_limit_factor )
    fprintf('\n GPU memory estimated / reco : %.2g MiB = %.2g GiB', gpu_mem_requ_per_reco/1024^2, gpu_mem_requ_per_reco/1024^3)    
    fprintf('\n GPU available memory / GiB:')
    fprintf(' %8.2g', mem_avail_gpu/1024^3)
    fprintf('\n number recos / GPU :')
    fprintf(' %7u', recos_per_gpu)
    fprintf('\n Parpool GPU index list:\n  ')
    fprintf(' %u', gpu_index_list)
    fprintf('\n Parpool GPU index list length: %u', numel(gpu_index_list))
    fprintf('\n')
end