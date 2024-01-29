classdef Transceiver
    %desired properties ['name', 'bandwidth', 'nFrequencySlots', 'GBaud', 'penalty', 'exampleCapacity']
    % 'name'            : [string] name string of a transceiver;
    % 'bandwidth'       : [1] signal bandwidth in GHz;
    % 'nFrequencySlots' : [1] number of occupied frequency slot;
    % 'GBaud'           : [1] baud rate of a transceiver;
    % 'Penalty'         : [1] SNR penalty in back-to-back configuration;
    % 'exampleCapacity' : [1] example transmission capaicty for this transceiver;
    %
    % Note:
    % 	- Show properties by closing function disp();
    properties
        name,
        bandwidth,          % in GHz, equals to GBaud by default.
        nFrequencySlots,	% minimumGrid
        GBaud,              % GBaud
        penalty,
        exampleCapacity,    % in Gbps, netBitRate,
    end
    methods
        
        
%         function str = disp(obj)
%             keyName = 'name';
%             if(numel(obj)>=2)
%                 str = sprintf('%s: [%s',keyName, obj(1).(keyName));
%                 for i = 2:numel(obj)
%                     str = sprintf('%s %s', str, obj(i).(keyName));
%                 end
%                 str = sprintf('%s]', str);
%             else
%                 str = sprintf('%s: %s',keyName, obj(1).(keyName));
%             end
%         end
        
        function obj = Transceiver(Vals)
            if nargin == 0
                return; % return empty obj.
            end;
            name            = Vals.('name');
            bandwidth       = Vals.('bandwidth');
            capacity        = Vals.('capacity');
            nFrequencySlots = Vals.('numberOfFrequencySlots');
            penalty         = Vals.('penalty');
            baudRate        = Vals.('baud_rate');
            assert(isnumeric(bandwidth)        ...
                &&isnumeric(nFrequencySlots)   ...
                &&isnumeric(capacity)          ...
                &&isnumeric(penalty)           ...
                &&isnumeric(baudRate)          ...
                , 'MATLAB:Transceiver:InputMustBeNumeric' ...
                , 'Value must be numeric');
            obj.name            = name;
            obj.bandwidth       = bandwidth;
            obj.nFrequencySlots = nFrequencySlots; % ParasTrans(idx_number_frequency_slots);
            obj.exampleCapacity = capacity;
            obj.penalty         = penalty;
            obj.GBaud           = baudRate;
        end
    end
end