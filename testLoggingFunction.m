function testDiary()

filename = '/mnt/tomoraid-LSDF/users/moosmann/testlog.txt';
delete(filename);
diary(filename);
try
    a = ones(10);
    fprintf('\n\nTEST LOGGING FUNCTION\n\nSize of a: %u x %u\n',size(a));
    a = a.*b;
    fprintf('\n\nSIze of %u %u\n',size(a));
    diary off
catch err
    report = getReport(err);
    fprintf('ERROR/EXCEPTION OCCURED:\n%s\n',report);
    diary off
end