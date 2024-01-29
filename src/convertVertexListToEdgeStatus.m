function isEdgeUsed = ...
    convertVertexListToEdgeStatus(...
    vertexList ...
    , edgeIndexMatrix)

    nEdge      = max(max(edgeIndexMatrix));
    nHop       = length(vertexList) - 1;
    isEdgeUsed = zeros(nEdge, 1);
    
    for iHop = 1 : nHop
        u = vertexList(iHop);
        v = vertexList(iHop + 1);

        isEdgeUsed(edgeIndexMatrix(u, v)) = true;
    end
end
