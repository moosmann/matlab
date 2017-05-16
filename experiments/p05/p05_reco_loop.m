function p05_reco_loop(nums, doreco, para)
% Script to loop over parameter sets related to paramters of scritp
% 'p05_reco'. Set parameters to loop over as elements of the structure
% array 'para' below. Fieldnames of 'para' MUST match the names of
% parameters in 'p05_reco'.
%
% Visual output ('visualOutput') and user interaction
% ('interactive_determination_of_rot_axis') are turned off by default if
% not set otherwise.
%
% Start loop by pushing 'F5', clicking on 'Run' in the Editor tab, or
% typing 'p05_reco_loop' in the Command Window.
%
% Written by Julian Moosmann. First version: 2017-02-15. Last: 2017-05-16

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nDATA SETS:')
for nn = 1:numel( para)
    [~, name] = fileparts( para(nn).scan_path );
    fprintf('\n%3u : %s', nn, name )
end

if ~isempty(nums)
    fprintf( '\n\nTO BE RECONSTRUCTED:')
    for nn = 1:numel( nums )
        num = nums(nn);
        external_parameter = para(num);
        [~, name] = fileparts( external_parameter.scan_path );
        fprintf('\n%3u : %s', num, name )                
    end
end

% Loop over parameter sets
if ~isempty( nums ) && doreco 
    fprintf( '\n\nSTART LOOPING \n')
    for nn = 1:numel( nums )
        num = nums(nn);
        
        external_parameter = para(num);
        [~, name] = fileparts( external_parameter.scan_path );
        fprintf('\nRECONSTRUCTION OF DATA SET NUMBER : %u : %s\n', num, name )
        
        p05_reco                
    end
    fprintf( '\nRECONSTRUCTION LOOP FINISHED')
end
fprintf('\n')
