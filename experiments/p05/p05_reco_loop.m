function p05_reco_loop(SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Script to loop over sets of paramters of reconstruction script 'p05_reco'
% as defined in global cell array 'PARAMETER_CELL' below.
%
% Visual output ('visual_output') and user interaction
% ('interactive_determination_of_rot_axis') are turned off by default if
% not set otherwise.
%
% Written by Julian Moosmann. First version: 2017-02-15. Last: 2017-06-04

global PARAMETER_CELL

try
    %% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 1
        SUBSETS = [];
    end
    if nargin < 2
        RUN_RECO = 0;
    end
    if nargin < 3
        PRINT_PARAMETERS = '';
    end
    
    %% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf( '\nDATA SETS:')
    for nn = 1:numel( PARAMETER_CELL)
        [~, name] = fileparts( PARAMETER_CELL{nn}.scan_path );
        fprintf('\n%3u : %s', nn, name )
    end
    
    if ~isempty(SUBSETS)
        fprintf( '\n\nTO BE RECONSTRUCTED:')
        for nn = 1:numel( SUBSETS )
            num = SUBSETS(nn);
            external_parameter = PARAMETER_CELL{num};
            [~, name] = fileparts( external_parameter.scan_path );
            if ~isempty(PRINT_PARAMETERS)
                if nn == 1
                    fprintf( '\n' )
                end
                if iscell(PRINT_PARAMETERS)
                    fprintf('%3u : %s\n', num, name )
                    for mm = 1:numel( PRINT_PARAMETERS )
                        fn = PRINT_PARAMETERS{mm};
                        fprintf( '        %s = ', fn)
                        %out = external_parameter.(fn);
                        out = eval(sprintf( 'external_parameter.%s', fn));
                        if isempty( out )
                            fprintf( '\n' )
                        else
                            if isstruct( out )
                                fprintf( '\n' )
                            end
                            disp( out )
                        end
                    end
                else
                    fprintf('%3u : %s', num, name )
                    fprintf( ', %s = ', PRINT_PARAMETERS)
                    %out = external_parameter.(PRINT_PARAMETERS);
                    fn = PRINT_PARAMETERS;
                    out = eval(sprintf( 'external_parameter.%s', fn));
                    if isempty( out )
                        fprintf( '\n' )
                    else
                        if isstruct( out )
                            fprintf( '\n' )
                        end
                        disp( out )
                    end
                end
            else
                fprintf('\n%3u : %s', num, name )
            end
        end
    end
    
    % Loop over parameter sets
    if ~isempty( SUBSETS ) && RUN_RECO
        fprintf( '\n\nSTART LOOPING \n')
        for nn = 1:numel( SUBSETS )
            num = SUBSETS(nn);
            
            external_parameter = PARAMETER_CELL{num};
            [~, name] = fileparts( external_parameter.scan_path );
            fprintf('\nRECONSTRUCTION OF DATA SET NUMBER %u : %s\n', num, name )
            
            p05_reco
        end
        fprintf( '\nRECONSTRUCTION LOOP FINISHED')
    end
    fprintf('\n')
    clear global
catch exception
    clear global
    rethrow( exception )
end
