function p05_reco_loop(SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Script to loop over sets of paramters of reconstruction script 'p05_reco'
% as defined in global cell array 'PARAMETER_CELL' below.
%
% Visual output ('visual_output') and user interaction
% ('interactive_mode.*') are turned off by default if
% not set otherwise.
%
% Written by Julian Moosmann. First version: 2017-02-15. Last: 2017-06-04

global PARAMETER_CELL

%try
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
    if ~iscell(PARAMETER_CELL{nn}.par.scan_path)
        [~, name] = fileparts( PARAMETER_CELL{nn}.par.scan_path );
        fprintf('\n%3u : %s', nn, name )
    else
        c = PARAMETER_CELL{nn}.par.scan_path;
        fprintf('\n%3u :', nn )
        for mm = 1:numel(c)
            [~, name] = fileparts(c{mm});
            fprintf(' %s',name )
        end
    end
end
if ~isempty(SUBSETS)
    fprintf( '\n\nTO BE RECONSTRUCTED:')
    for nn = 1:numel( SUBSETS )
        num = SUBSETS(nn);
        external_parameter = PARAMETER_CELL{num};
        c = external_parameter.par.scan_path;
        [~, name] = fileparts( external_parameter.par.scan_path );
        if ~isempty(PRINT_PARAMETERS)
            if nn == 1
                fprintf( '\n' )
            end
            if iscell(name)
                kke = numel(name);
            else
                kke = 1;
            end
            for kk = 1:kke
                if iscell(name)
                    name_kk = name{kk};
                else
                    name_kk = name;
                end
                if iscell(PRINT_PARAMETERS)
                    
                    fprintf('%3u : %s\n', num, name_kk )
                    for mm = 1:numel( PRINT_PARAMETERS )
                        fn = PRINT_PARAMETERS{mm};
                        fprintf( '        %s = ', fn)
                        out = eval(sprintf( 'external_parameter.%s', fn));
                        if isempty( out )
                            fprintf( '\n' )
                        else
                            if isstruct( out )
                                fprintf( '\n' )
                            end
                            fprintf('%s\n',num2str(out))
                        end
                    end
                else
                    fprintf('%3u : %s', num, name_kk )
                    fprintf( ', %s = ', PRINT_PARAMETERS)
                    fn = PRINT_PARAMETERS;
                    out = eval(sprintf( 'external_parameter.%s', fn));
                    if isempty( out )
                        fprintf( '\n' )
                    else
                        if isstruct( out )
                            fprintf( '\n' )
                        end
                        fprintf('%s',num2str(out))
                    end
                end
            end
        else
            if ~iscell(name)
                fprintf('\n%3u : %s', num, name )
            else
                fprintf('\n%3u :',num)
                for mm = 1:numel(c)
                    [~, name] = fileparts(c{mm});
                    fprintf(' %s',name )
                end
            end
        end
    end
end

% Loop over parameter sets
if ~isempty( SUBSETS ) && RUN_RECO == 1
    fprintf( '\n\nSTART LOOPING \n')
    for nn = 1:numel( SUBSETS )
        num = SUBSETS(nn);
        external_parameter = PARAMETER_CELL{num};
        c = external_parameter.par.scan_path;
        [~, name] = fileparts( external_parameter.par.scan_path );
        if ~iscell(name)
            fprintf('\nRECONSTRUCTION OF DATA SET NUMBER %u : %s\n', num, name )
        else
            fprintf('\nRECONSTRUCTION OF DATA SET NUMBER %u :', num )
            for mm = 1:numel(c)
                [~, name] = fileparts(c{mm});
                fprintf(' %s',name )
            end
            fprintf('\n')
        end

        p05_reco(external_parameter)
    end
    fprintf( '\nRECONSTRUCTION LOOP FINISHED')
end
fprintf('\n')
clear global
%     catch exception
%         clear global
%         rethrow( exception )
% end
