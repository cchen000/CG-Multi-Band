function costMatrix = fixCostMatrix(distanceMatrix)
    %(i,j) : 0 -> inf(disconnected)
    % 
    % Input:
    %       distanceMatrix      distance between nodes i and j;
    % Output:
    %       costMatrix          network routing cost between nodes i and j;
    
    nNode       = size(distanceMatrix, 1);
    costMatrix  = zeros(nNode, nNode);
    
    for iNode = 1 : nNode
        for jNode = 1 : nNode
            if iNode ~= jNode
                if distanceMatrix(iNode, jNode) == 0
                    costMatrix(iNode, jNode) = inf;
                else
                    costMatrix(iNode, jNode) = distanceMatrix(iNode, jNode);
                end
            else
                costMatrix(iNode, jNode) = distanceMatrix(iNode, jNode);
            end
        end
    end
end