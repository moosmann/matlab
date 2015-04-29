%% stage 11, main emby
TomoTaken = 1:14;
VolCor = [3 5 10 13];
TomoUsed = TomoTaken;
TomoUsed(VolCor) = [];
deadtime = 10;
tomotime = (50+1200+50+60+60)*0.015; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 11, main emby, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
fprintf('Corrupted volumes :')
fprintf('%5u ',VolCor)
fprintf('\n')
fprintf('Volumes used:  ')
fprintf('%5u ',TomoUsed)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n\n')

%% stage 14
TomoTaken = 1:10;
VolCor = [1 3 9];
TomoUsed = TomoTaken;
TomoUsed(VolCor) = [];
deadtime = 10;
tomotime = (50+1200+50+60+60)*0.025; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 14,SupVideo 5,tomo time: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
fprintf('Corrupted volumes :')
fprintf('%5u ',VolCor)
fprintf('\n')
fprintf('Volumes used:  ')
fprintf('%5u ',TomoUsed)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n\n')


%% stage 12.5 
TomoTaken = 1:8;
TomoUsed = 1:7;
deadtime = 10; % min
tomotime = (50+1200+50+60+60)*0.05; % s
tempRes = tomotime + deadtime*60;
fprintf('STAGE 12.5, suppl movie 4b, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')

%% October data, Stage 11
TomoTaken = (0:10)+1;
TomoUsed = TomoTaken;
deadtime = 8; % min
tomotime = 50; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 11, October data, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')

%% October data, Stage 11
TomoTaken = (0:8)+1;
TomoUsed = TomoTaken;
deadtime = 10; % min
tomotime = 50; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 11, October data, missing frames (externally saved), time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')

%% stage 16
TomoTaken = 1:5;
TomoUsed = TomoTaken;
deadtime = 5; % min
tomotime = (50+1200+50+60+60)*0.05; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 16, SupplVideo 9, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')

%% stage 18
TomoTaken = 1:8;
TomoUsed = 1:7;
deadtime = 5; % min
tomotime = (50+1200+50+60+60)*0.075; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 18, SupplVideo 10, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')

%% stage 11.5
TomoTaken = 1:8;
TomoUsed = 1:8;
deadtime = 2; % min
tomotime = (50+1200+50+60+60)*0.05; % s
tempRes = tomotime + deadtime*60;
fprintf('\nSTAGE 11.5, unreleased, time for tomogram: %g, temporal resolution: %g\n',tomotime,tempRes/60);
TomoTimePts = (tomotime/2 + (deadtime*60+tomotime)*(TomoTaken-1))/60;
% Print
fprintf('Volumes taken: ')
fprintf('%5u ',TomoTaken)
fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.1f ',TomoTimePts)
fprintf('\n')
% fprintf('Amira counting:')
% fprintf('%5u ',0:numel(TomoUsed)-1)
% fprintf('\n')
fprintf('  Time in min: ')
fprintf('%5.0f ',TomoTimePts(TomoUsed))
fprintf('\n')
