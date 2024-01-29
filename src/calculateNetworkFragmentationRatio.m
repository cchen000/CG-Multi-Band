function FragRatio = calculateNetworkFragmentationRatio(isEdgeUseSlot, bandwidthTypesSet, FormatOptions)
%calculateNetworkFragmentationRatio Calcualte network basic informations.
% 
% Derived from : D:\Git\test_frequency\code\CaluclateEntropyMetrics
% Derived Date : 15, Jun. 2022.
% v2: add simple and detailed output printouts. date: 15, Jun. 2022. by cchen
% --------
% Derived from D:\Git\Elastic_Service_Journal\code\CalculateNetworkInfos\CalculateNetworkInfos.m
% Date: 5th, Mar. 2023.
FormatOptions.SpectrumFragmentation = 'Shannon';
[NoEdges, NUMBER_OF_SLOTS] = size(isEdgeUseSlot);

% ==============================
% - Display spectrum fragmentaion ratio;
% ==============================
H = 0;
switch(FormatOptions.SpectrumFragmentation)
    case 'Shannon' % Shannon entropy
        for lth = 1:NoEdges
            H = H + calculateFragmentaionShannon(logical(isEdgeUseSlot(lth,1:NUMBER_OF_SLOTS)));
        end
    case 'ABP' % Access blocking probability
        for lth = 1:NoEdges
            H = H + calculateFragmentaionABP(logical(isEdgeUseSlot(lth,1:NUMBER_OF_SLOTS)), bandwidthTypesSet);
        end
    otherwise
        error('none specifie');
end
if(isnan(H))
    error('NaN error');
end
fprintf('Spectrum fragmentation ratio = %g\n',H);

% ==============================
% - Return parameter
% ==============================
if nargout==1
    FragRatio = H;
end

