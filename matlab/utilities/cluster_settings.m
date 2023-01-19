c = parcluster();
c.AdditionalProperties.Username = 'moosmanj';
c.AdditionalProperties.QueueName = 'maxcpu';
c.AdditionalProperties.WallTime = '00:05:00';
c.AdditionalProperties.EmailAddress = 'julian.moosmann@hereon.de';
%c.AdditionalProperties.UseIdentityFile =true;
%c.AdditionalProperties.IdentityFile = '/home/moosmanj/.ssh/id_rsa';
c.AdditionalProperties.UseIdentityFile = false;
c.AdditionalProperties.IdentityFile = '';
c.AdditionalProperties.IdentityFileHasPassphrase = false;
c.AdditionalProperties.ClusterHost = hostname;
c.ClusterMatlabRoot = matlabroot;
c.AdditionalProperties.EnableDebug = true;
c.saveProfile;


cdhome
job = c.batch(@pwd, 1, {}, 'CurrentFolder','.');



