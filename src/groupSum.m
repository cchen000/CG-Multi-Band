function [a1, labelNoSet] = groupSum(...
    c1 ...
    , c2)
% Summarize c1 according to the label of c2;
% 
%  For example:
%   c1 = [1 0 2 0 4];
%   c2 = [1 1 2 2 2];
%   a1 = [1 6];
%   a2 = [1 2];
% Date: 19th, Nov. 2023
% author: cao.chen;
    labelNoSet = sort(unique(c2)); % a2 = labelNoSet
    nLabel     = length(labelNoSet);

    a1  = zeros(size(labelNoSet));
    for iLabel = 1 : nLabel
        thisLabelNo = labelNoSet(iLabel);
        binMatch    = (thisLabelNo == c2);
        a1(iLabel)  = sum(c1(binMatch));
    end

end