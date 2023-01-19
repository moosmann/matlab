function write_hdf(filename, arrray, datasetname)

if nargin < 3
    datasetname = ['/' inputname(2)];
end

a_class = class(array);
a_size = size(array);

h5create(filename,datasetname,a_size,'Datatype',a_class)
h5write(filename,datasetname,arrray)
