
%% Example case - DT9
clear all;
ImportLib();

% ==============================
% Read simulation configurations
% ==============================
% Default setting 
% ::details are either given in .yaml file;
configFile                          = 'setup/SimulationSetup_DT9.yaml';
SimulationSetting                   = ReadYaml(configFile); 

[LightPathSet, isEdgeUseSlot, info] = main(SimulationSetting);

fprintf("Running time = %g seconds\n",info('Time'));

pause(0.1);
figure(87),
clf;
imagesc(logical(isEdgeUseSlot), 'CDataMapping', 'scaled');colorbar
title('Spectrum Usage Graph');
xlabel('Slot');
ylabel('Edge');
hold off;
