
clear all
%% Fast switch
par.quick_switch = 1;
par.a = 1;
par.b = 2;
tomo.a = 1;
par_quick_switch.who = who();
par_quick_switch.num_structs = numel( par_quick_switch.who );
for n = 1:par_quick_switch.num_structs
    sn = par_quick_switch.who{n};
    par_quick_switch.structs.(sn) = eval( sn );
end



%% Parameters
par.a = 11;
par.b = 22;
par.c = 30;
tomo.a = 11;
tomo.b = 22;
tomo.c = 33;
% %%
if par.quick_switch
    
    for ss = 1:par_quick_switch.num_structs
        sn = par_quick_switch.who{ss};        
        tmp = par_quick_switch.structs.(sn);
        
        fn = fieldnames( tmp);
        for nn = 1:numel( fn )            
            field = fn{nn};            
            eval_string = sprintf( '%s.(''%s'') = par_quick_switch.structs.(''%s'').(''%s'');', sn, field, sn, field );            
            eval( eval_string );
        end
        
        
    end
end

disp( 'FINAL PARAMETERS' )
par
tomo



