function PrePath(String)
% Check hostname and modify the input path 'String' according the machines
% mount of the LSDF storage

[~, hostname] = unix('echo $HOSTNAME');
hostname(end-1:end) = [];
switch hostname
    case 'anka-nc-gpu'
        String = ['/mnt/gpuNodeStorage/scratch/moosmann/' String];
        %ParentPath = ['/mnt/LSDF/tomo/' ParentPath];
    case 'anka-nc-cpu'
        String = ['/mnt/LSDF/tomo/' String];
    otherwise
        String = ['/mnt/tomoraid-LSDF/tomo/' String];
end

assignin('caller',inputname(1),String)
