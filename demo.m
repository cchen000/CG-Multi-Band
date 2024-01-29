%% Example for changing parameters
%% 
% 
%  Example for changing optical bands

SimulationSetting.number_of_bands               = 3;
SimulationSetting.widebandset                   = 'UBand, CBand, LBand ';
SimulationSetting.wideband_file                 = 'networks\raw\DT9_ULCBands.txt';

%  SimulationSetting.number_of_bands     = 1;
%  SimulationSetting.widebandset         = 'ULCBand ';
%  SimulationSetting.wideband_file       = 'networks\raw\DT9_ULCBand.txt';

%%
%  Example for changing TRx
SimulationSetting.number_of_transponder = 1;           % currently supporting only 1 baud-rate ;
SimulationSetting.transponder_set       = 'Trans50G'; % 500 GBaud TRx.
SimulationSetting.frequency_slots       = 50;         % How network spectrum resource is divided
SimulationSetting.number_of_slots       = 15000/50;

%%
%  Example for changing transmission modes
% 
SimulationSetting.number_of_transmission_modes  = 8;
SimulationSetting.transmission_mode_set         = 'PM_BPSK, PM_QPSK, PM_8QAM, PM_16QAM, PM_32QAM, PM_64QAM, PM_128QAM, PM_256QAM';
