weblink_url = 'https://photon-science.desy.de/facilities/petra_iii/machine/parameters/index_eng.html';
weblink_name = 'photon-science.desy.de/facilities/petra_iii/machine/parameters';
weblink = sprintf('<a href = "%s">%s</a>', weblink_url, weblink_name);

energy = 6e9; %eV
r = 2304; %m
num_buckets = 3840;
modes = {'normal', 'multi', 'timing'};
num_bunches = [960, 480, 40]; % multi bunch, multi bunch, timing mode
bucket_filling = num_buckets ./ num_bunches;
bunch_separation = [8, 16, 192]*1e-9; %s
rev_time = 7.685e-6; %s
rev_freq = 1 / rev_time; % Hz
rev_time_c = r / SpeedOfLight;

electron_beam_current = 100e-3; %A
electron_beam_charge = 769e-9; %C
electron_beam_lifetime = [13, 13 1] * 3600; %s
num_electrons = 4.8e12;
num_electrons_from_beam_charge = electron_beam_charge / ElementaryCharge;

bunch_length = 0.0132; %m

%% Experiment

frame_rate = 1e3; % Hz
exp_time = 1/ frame_rate;
bunches_per_exptime = exp_time ./ bunch_separation;

%% Print Info

fprintf( 'PETRA III machine parameters' )

fprintf( '\nweb: %s', weblink )
disp( weblink )

fprintf( 'energy: %0.f keV', energy / 1e9 )
fprintf( '\ncircumference: %f m', r )
fprintf( '\nnumber of buckets: %u', num_buckets )
fprintf( '\nrevolution time: %.4f µs', rev_time * 1e6 )
fprintf( '\nrevolution time at light speed: %.4f µs', rev_time_c * 1e6 )
fprintf( '\nrevolution frequency: %.1f kHz', rev_freq/1000 )

fprintf( '\nelectron beam current: %.0f mA', electron_beam_current*1000 )
fprintf( '\nelectron beam charge: %.0f nC', electron_beam_charge*1e9 )
fprintf( '\nelectron beam current * ref_time: %0.f nC', 1e9*electron_beam_current * rev_time )

fprintf( '\nnumber of electrons:      %10g', num_electrons )
fprintf( '\nelectron beam charge / e: %10g', num_electrons_from_beam_charge )

fprintf( '\nbunch length (rms): %.1f mm or %.1f ps', bunch_length * 1000, bunch_length / SpeedOfLight * 1e12 )

fprintf( '\nExperimental parameters' )
fprintf( '\nframe rate: %.3f kHz', frame_rate / 1000 )
fprintf( '\nexposure time: %.2f ms', exp_time * 1000 )

fprintf( '\nMACHINE' )
fprintf( '\n%30s: ', 'modes' )
fprintf( '%10s', modes{:})

fprintf( '\n%30s: ','bucket filling pattern')
fprintf( '%10u', bucket_filling )

fprintf( '\n%30s: ','number of bunches')
fprintf( '%10u', num_bunches )

fprintf( '\n%30s: ','bunch separation / ns')
fprintf( '%10g', bunch_separation * 1e9 )

fprintf( '\n%30s: ', 'revolution time / #bunches' )
fprintf( '%10.1f', rev_time * 1e9 ./ num_bunches )

fprintf( '\n%30s: ', 'electron beam lifetime / h' )
fprintf( '%10.0f', electron_beam_lifetime / 3600)

fprintf( '\n%30s: ', '1% lifetime / min' )
fprintf( '%10.0f', 0.01 * electron_beam_lifetime / 60)

fprintf( '\n%30s: ', 'time for 1% drop / min' )
fprintf( '%10.0f', -log(0.99) * electron_beam_lifetime / 60 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nEXPERIMENT' )
fprintf( '\n%30s: ','bunches / exposure time' )
fprintf( '%10.0f', bunches_per_exptime )


fprintf( '\n' )