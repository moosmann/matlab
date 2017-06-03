function p05_reco_loop(nums, doreco, print_field, par)
% Script to loop over parameter sets related to paramters of script
% 'p05_reco'. Set parameters to loop over as elements of the structure
% array 'par' below. Fieldnames of 'par' MUST match the names of
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

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CELL ARRAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell( par )
    
    fprintf( '\nDATA SETS:')
    for nn = 1:numel( par)
        [~, name] = fileparts( par{nn}.scan_path );
        fprintf('\n%3u : %s', nn, name )
    end
    
    if ~isempty(nums)
        fprintf( '\n\nTO BE RECONSTRUCTED:')
        for nn = 1:numel( nums )
            num = nums(nn);
            external_parameter = par{num};
            [~, name] = fileparts( external_parameter.scan_path );
            if ~isempty(print_field)
                if nn == 1
                    fprintf( '\n' )
                end
                if iscell(print_field)
                    fprintf('%3u : %s\n', num, name )
                    for mm = 1:numel( print_field )
                        fn = print_field{mm};
                        fprintf( '        %s = ', fn)
                        out = external_parameter.(fn);
                        if isempty( out )
                            fprintf( '\n' )
                        else
                            disp( out )
                        end
                    end
                else
                    fprintf('%3u : %s', num, name )
                    fprintf( ', %s = ', print_field)
                    out = external_parameter.(print_field);
                    if isempty( out )
                        fprintf( '\n' )
                    else
                        disp( out )
                    end
                end
            else
                fprintf('\n%3u : %s', num, name )
            end
        end
    end
    
    % Loop over parameter sets
    if ~isempty( nums ) && doreco
        fprintf( '\n\nSTART LOOPING \n')
        for nn = 1:numel( nums )
            num = nums(nn);
            
            external_parameter = par{num};
            [~, name] = fileparts( external_parameter.scan_path );
            fprintf('\nRECONSTRUCTION OF DATA SET NUMBER %u : %s\n', num, name )
            
            p05_reco
        end
        fprintf( '\nRECONSTRUCTION LOOP FINISHED')
    end
    fprintf('\n')

    
%% STRUCT ARRAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif isstruct( par )
    
    fprintf( '\nDATA SETS:')
    for nn = 1:numel( par)
        [~, name] = fileparts( par(nn).scan_path );
        fprintf('\n%3u : %s', nn, name )
    end
    
    if ~isempty(nums)
        fprintf( '\n\nTO BE RECONSTRUCTED:')
        for nn = 1:numel( nums )
            num = nums(nn);
            external_parameter = par(num);
            [~, name] = fileparts( external_parameter.scan_path );
            if ~isempty(print_field)
                if nn == 1
                    fprintf( '\n' )
                end
                if iscell(print_field)
                    fprintf('%3u : %s\n', num, name )
                    for mm = 1:numel( print_field )
                        fn = print_field{mm};
                        fprintf( '        %s = ', fn)
                        out = external_parameter.(fn);
                        if isempty( out )
                            fprintf( '\n' )
                        else
                            disp( out )
                        end
                    end
                else
                    fprintf('%3u : %s', num, name )
                    fprintf( ', %s = ', print_field)
                    disp( external_parameter.(print_field) )
                end
            else
                fprintf('\n%3u : %s', num, name )
            end
        end
    end
    
    % Loop over parameter sets
    if ~isempty( nums ) && doreco
        fprintf( '\n\nSTART LOOPING \n')
        for nn = 1:numel( nums )
            num = nums(nn);
            
            external_parameter = par(num);
            [~, name] = fileparts( external_parameter.scan_path );
            fprintf('\nRECONSTRUCTION OF DATA SET NUMBER %u : %s\n', num, name )
            
            p05_reco
        end
        fprintf( '\nRECONSTRUCTION LOOP FINISHED')
    end
    fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    fprintf( '\n INVALID PARAMETER STRUCT! \n')
end
