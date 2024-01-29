function [configNoThisBand, repConfigNo, repTimes] = ...
    sortConfig(...
    configNo ...
    , repeatConfigurations ...
    , sortMetric...
    , direction...
    )
% 
% Description
%       Repeat configuration by configuration times
% 
% For example:
%       ConfigNo    = [1 2 3 4 5 6];
%       repeatConfigurations = [1 0 2 0 0 3];
%       sort_metric = [1 0 2 0 0 3]; % repeatConfigurations
%     ::repConfigNo = [6 3 1];
%     ::repTimes    = [3 2 1];
%     =>configNoThisBand = [6 6 6 3 3 1];
%     Or,
%       ConfigNo    = [1 2 3 4 5 6];
%       repeatConfigurations = [1 0 2 0 0 3]; 
%       sort_metric = [20 1 10 100 100 1000];
%     ::repConfigNo = [6 1 3];
%     ::repTimes    = [3 1 2];
%     =>configNoThisBand = [6 6 6 1 3 3];
%     Or,
%       ConfigNo    = [1 2 3 4 5 6];
%       repeatConfigurations = [1 0 2 0 0 3]; 
%       sort_metric = [20 1 1000 100 100 10];
%     ::repConfigNo = [3 1 6];
%     ::repTimes    = [2 1 3];
%     =>configNoThisBand = [3 3 1 6 6 6];
%     Or,
%       ConfigNo    = [1 6 5 3 2 4];
%       repeatConfigurations = [1 3 0 2 0 0]; 
%       sort_metric = [20 1 1000 100 100 10];
%     ::repConfigNo = [3 1 6];
%     ::repTimes    = [2 1 3];
%     =>configNoThisBand = [3 3 1 6 6 6];
% 
%  Date: 20th, Nov. 2023.
%  Author: cao chen;

    if nargin == 3
        direction = 'descend';
    else
        % deal with four inputs
    end

    % find active configuration;
    isActive_ofconfig = logical(repeatConfigurations);

    % sort rows;
    [val] = sortrows(...
        [sortMetric(:), isActive_ofconfig(:), configNo(:), repeatConfigurations(:)], ...
        1, direction);

    % find previous activating status after sorting;
    isStillActive_ofconfig = logical(val(:, 2));

    % get parameters after sorted;
    repTimes    = val(isStillActive_ofconfig, 4)';
    repConfigNo = val(isStillActive_ofconfig, 3)';

    % repeat parameters;
    configNoThisBand = repelem(repConfigNo, repTimes);
end