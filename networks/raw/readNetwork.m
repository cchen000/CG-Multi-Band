function [varargout] = readNetwork(fileName, format)
%
% Description:
%       Read network topology profiles
%
% Input:
% 	FileName : network topology information;
%	format : personlalized network label.
% Output:
%   [nNodes, nEdges, distanceMatrix, longitudeDegreeArray, latitudeDegreeArray]
%
% See also: ImportNetworkFromFile.m, and notes.txt
%
% Date: 4th, Dec. 2023.
% Author: cao.chen
% ==============================
% Read the following information:
% :: basic
% :: edge
% :: node

if nargin==0
    % Example 1;
    fileName = 'DT14.txt';
    %     'DT9.txt';
    format =struct(...
        'Direction', 'directed', ...
        'LinkRepresentation', 'edgeAdjacent', ...
        'NodeCoordinate', true);
    
%     % Example 2;
%     fileName = 'format2_DT14.txt';
%     %     'format2_DT9.txt';
%     %     'format2_EX4.txt';
%     format =struct(...
%         'Direction', 'directed', ...
%         'LinkRepresentation', 'distanceMatrix', ...
%         'NodeCoordinate', false);
    
    
%     % Example 3;
%     fileName = 'format3_DT9.txt';
%     %     'format2_DT9.txt';
%     %     'format2_EX4.txt';
%     format =struct(...
%         'Direction', 'undirected', ...
%         'LinkRepresentation', 'distanceMatrix', ...
%         'NodeCoordinate', false);
end

% open file
fid = fopen(fileName,'r+');
assert(fid>=1, 'Wrong Reading File');

% ==============================
% Read basic information
% ==============================
% :: nNodes;
% :: nEdges;
[A,COUNT] = fscanf(fid,'%d %d',[1,2]);
nNode   = A(1); 
nEdges  = A(2);

% ==============================
% Read edge information;
% ==============================
distanceMatrix = zeros(nNode,nNode);
if(strcmp(format.('LinkRepresentation'),'edgeAdjacent'))
    % Type 1
    FormatSpec = repmat('%d %d %d',[1,1]);
    for i=1:nEdges
        [A,COUNT] = fscanf(fid,FormatSpec,[1,3]);
        distanceMatrix(A(1),A(2) ) = A(3);
    end
elseif(strcmp(format.('LinkRepresentation'),'distanceMatrix'))
    % Type 2
    FormatSpec = repmat('%d ',[1,nNode]);
    for i=1:nNode
        [A,COUNT] = fscanf(fid,FormatSpec,[1,nNode]);
        distanceMatrix(i,1:COUNT) = A;
    end
else
    error('not specified');
end

% Check edge direction;
if(strcmp(format.('Direction'),'directed'))
    ;
elseif(strcmp(format.('Direction'),'undirected'))
    for mRow = 1:nNode
        for nCol = 1:nNode
            % Check Duv==Dvu if Duv or Duv are not 0.
            assert(distanceMatrix(mRow,nCol)*distanceMatrix(nCol,mRow)==0 ...
                ,'not symmetric');
        end 
    end
    % Assign bidirectional possiblity for undircted graph;
    distanceMatrix      = distanceMatrix + distanceMatrix';
end

% ==============================
% Read node information
% ==============================
latitudeDegreeArray = zeros(nNode,1);
longitudeDegreeArray= zeros(nNode,1);
if(format.('NodeCoordinate')==true)
    FormatSpec          = repmat('%d %g %g',[1,1]);
    for i=1:nNode
        [A,COUNT] = fscanf(fid,FormatSpec,[1,3]);
        longitudeDegreeArray(A(1))  = A(2);% east and west
        latitudeDegreeArray(A(1))   = A(3);% south and north
    end
elseif(format.('NodeCoordinate')==false)
    ; % do nothing about node coordinates.
end
fclose(fid);

% ==============================
% Set output.
% ==============================
        nOutputs = nargout;
        varargout = cell(1,nOutputs);
        
        for k = 1:nOutputs
            switch(k)
                case 1
                    varargout{1} = nNode;
                case 2
                    varargout{2} = nEdges;
                case 3
                    varargout{3} = distanceMatrix;
                case 4
                    varargout{4} = longitudeDegreeArray;
                case 5
                    varargout{5} = latitudeDegreeArray;
                otherwise
                    error('undefined');
            end
        end
end
