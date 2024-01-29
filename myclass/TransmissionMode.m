classdef TransmissionMode
% including name, code rate, capacity, min. required SNRs, and max. reach.
% 
%     'name'              : ['string'] e.g., PM_BPSK
%     'exampleClientRate' : [1] example client bit rate;
%     'SNRdB'             : [1] minimum required SNR in decibel
%     'SNR'               : [1] minimum required SNR
%     'maxReach'          : [1] default maximum transmission reach
%
%     'spectralEfficiency': [1] spectral efficiency.
%
%
% Close display by comment disp();
    properties
        name    = '';
        SNRdB   = [];   % [dB]
        SNR     = [];
        
        exampleClientRate  = []; % [Gbps]
        spectralEfficiency = []; % [bps/Hz]
        maxReach           =  0; % [span]
    end
    methods
        
        
        function str = disp(obj)
            str = dispName(obj);
        end
        function obj = TransmissionMode(Vals)
            if nargin == 0
                return;
            end 
            name     = Vals.('name');
            SNRdB    = Vals.('SNRdB');
            dataRate = Vals.('dataRate');
            maxReach = Vals.('reach');

            defaultBaud = Vals.('defaultBaud');
            
            assert(isnumeric(SNRdB) ...
                && isnumeric(dataRate) ...
                && isnumeric(maxReach)...
                , 'MATLAB:TransmissionMode:InputMustBeNumeric' ...
                , 'Value must be numeric');
            
            obj.name    = name;
            obj.SNRdB   = SNRdB;
            obj.SNR     = dB2lin( SNRdB );
            
            obj.exampleClientRate = dataRate;
            obj.maxReach = maxReach;
            
            obj.spectralEfficiency = ...
                round(dataRate/defaultBaud,1);
        end
    end
end