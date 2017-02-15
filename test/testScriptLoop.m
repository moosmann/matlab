function testScriptLoop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over
para(1).scan_path = 'path1';
para(1).para1 = 22;

para(2).scan_path = 'path2';
para(2).scan_para2 = 23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop over all parameter sets
for nn = 1:numel( para )
    
    external_parameter = para(nn);    
    testScriptToLoop
end