function ImportLib()

% add path;
addpath(genpath('src\'));
addpath(genpath('myclass\'));

% Check function;
if(     ~exist('ReadYaml.m') ||...      % Missing ReadingYAML.m
        ~exist('OpticalBand.m') || ... % Missing optical band;
        ~exist('LightPath.m') || ... % Missing lightpath;
        ~exist('Demand.m') || ... % Missing demand;
        ~exist('TransmissionMode.m') || ... % Missing the modulation format;
        ~exist('Transceiver.m')...    % Missing the transponders;
        )
    error('main: Files missing');
end

end