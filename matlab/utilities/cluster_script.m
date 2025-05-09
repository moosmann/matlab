function out = cluster_script()

pause(2);
hn = getenv('HOSTNAME');
p = '/home/moosmanj/slurm';
t = datetime('now','Format','yyyy-MM-dd_hh-mm-ss-SSS');
fn = sprintf('%s/%s_%s.txt',p,hn,t);
%fprintf('\n%s', fn)
fid = fopen(fn,"w");
fwrite(fid,fn);
%fwrite(fid,'\n');
fwrite(fid,pwd);
fclose(fid);

out = 1;