function ratio = calculateNetworkOccupationRate(isEdgeUseSlot)
%ratio = total used slots  / total slots.
    [nEdges, nSlots] = size(isEdgeUseSlot);
    ratio            = sum(sum(isEdgeUseSlot)) / nEdges / nSlots;
    fprintf('Network Occuptation Rate = %g \n', ratio);
end