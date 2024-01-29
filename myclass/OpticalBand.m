classdef OpticalBand < matlab.mixin.Copyable % allow copyable property. 
% OpticalBand ('WideOpticalBand')
%      'name'               : ['string'] e.g., LBand, CBand;
%      'spectrum'           : [1] bandwidth of an optical band [GHz];
%      'nFrequencySlots'    : [1] number of slots;
%      'limitedTransceiver' : [1] maximum allowed transceiver
%      'SNRpenaltydB'       : [1] SNR margin [dB];
%      'SNRpenalty'         : [1] SNR margin;
%      'sbegin'             : slot beginning;
%      'send'               : slot ending;
%      'centerFrequency'    : relative center frequency of a band [GHz].
%
% Close display by commenting function disp();
    properties
        name                  = '';
        spectrum              =  0;
        nFrequencySlots       =  0;
        nAvailabelTransceiver =  0;
        SNRpenaltydB          = -1;
        SNRpenalty            = -1;
        sbegin                =  0;
        send                  =  0;
        centerFrequency       =  0; 
    end
    methods
        
        function str = disp(obj)
            KeyName = 'name';
            if(numel(obj)>=2)
                str = sprintf('%s: [%s',KeyName, obj(1).(KeyName));
                for i = 2:numel(obj)
                    str = sprintf('%s %s', str, obj(i).(KeyName));
                end
                str = sprintf('%s]', str);
            else
                str = sprintf('%s: %s',KeyName, obj(1).(KeyName));
            end
        end
        
        
        function obj = OpticalBand(Vals)
            if(nargin==0)
                return;
            end
            bandName             = Vals.('name');
            spectrumResource     = Vals.('spectrum_resource');
            nAvailabelTransceiver= Vals.('max_transceivers');
            SNRpenaltydB         = Vals.('SNRpenaltydB');
            centerFrequency      = Vals.('center_frequency');
            frequencySlots       = Vals.('frequency_slots');
            
            assert(isnumeric(spectrumResource), 'Value must be numeric');
            assert(isnumeric(nAvailabelTransceiver), 'Value must be numeric');
            assert(isnumeric(SNRpenaltydB), 'Value must be numeric');
            assert(isnumeric(centerFrequency), 'Value must be numeric');
            assert(isnumeric(frequencySlots), 'Value must be numeric');

            obj.name                  = bandName;
            obj.spectrum              = spectrumResource;
            obj.nFrequencySlots       = spectrumResource / frequencySlots;
            obj.nAvailabelTransceiver = nAvailabelTransceiver;
            obj.SNRpenaltydB          = SNRpenaltydB;
            obj.SNRpenalty            = dB2lin(SNRpenaltydB);
            obj.centerFrequency       = centerFrequency;
        end
    end
end