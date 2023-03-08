c = parcluster();
c.AdditionalProperties.AdditionalSubmitArgs = '--nodes=5';
c.AdditionalProperties.Username = 'moosmanj';
c.AdditionalProperties.QueueName = 'maxcpu';
c.AdditionalProperties.WallTime = '00:05:00';
%c.AdditionalProperties.EmailAddress = 'julian.moosmann@hereon.de';
c.AdditionalProperties.UseIdentityFile =true;
c.AdditionalProperties.IdentityFile = '/home/moosmanj/.ssh/id_rsa';
%c.AdditionalProperties.UseIdentityFile = false;
%c.AdditionalProperties.IdentityFile = '';
c.AdditionalProperties.IdentityFileHasPassphrase = false;
c.AdditionalProperties.ClusterHost = getenv('HOSTNAME');
%c.ClusterMatlabRoot = matlabroot;
%c.AdditionalProperties.EnableDebug = true;c.A
c.saveProfile;

%job = c.batch(@pwd, 1, {}, 'CurrentFolder','.');
p = '/home/moosmanj/slurm/';
job = c.batch(@cluster_script, 1, {}, 'CurrentFolder',p);



