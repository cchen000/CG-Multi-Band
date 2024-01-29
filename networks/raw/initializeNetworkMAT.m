function initializeNetworkMAT(importFile, exportMat)
% Description:
%       create .mat format file for networks.
% 
% [importFile] topo file with edge information and node information;
% [exportMat]  name of .mat files to be saved.
% 
% 
% Date: 6th, Dec., 2023
% Author: cao chen
% ==============================
% load file;
% ==============================
if(nargin==0)
    importFile = 'EX4.txt';
    exportMat  = 'EX4.mat';
end
format = struct(...
    'Direction', 'directed', ...
    'LinkRepresentation', 'edgeAdjacent', ...
    'NodeCoordinate', true);

[nNodes, nEdges, distanceSpanMatrix, longitudeDegreeArray, latitudeDegreeArray]  = readNetwork(importFile, format);

% ==============================
% Calculate
% ==============================

% Euclidian distance;
absEuclidianMatrix = zeros(nNodes,nNodes);
for i = 1:(nNodes*(nNodes-1))
    [s,d] = ind2sub([nNodes,nNodes],i);
    
    lat1 = latitudeDegreeArray(s)*pi/180;
    lat2 = latitudeDegreeArray(d)*pi/180;
    lon1 = longitudeDegreeArray(s)*pi/180;
    lon2 = longitudeDegreeArray(d)*pi/180;
    Point1 = struct('latitude',lat1, ...
        'longitude', lon1);
    Point2 = struct('latitude',lat2, ...
        'longitude', lon2);
    absEuclidianMatrix(s,d) = fun_harversine(Point1,Point2);
    absEuclidianMatrix(d,s) = absEuclidianMatrix(s,d);
end

% Fiber distance;
absFiberMatrix   = zeros(nNodes,nNodes);
for i = 1:(nNodes*(nNodes-1))
    [s,d] = ind2sub([nNodes,nNodes],i);
    if(distanceSpanMatrix(s,d)>=1e-4)
        absFiberMatrix(s,d) =  estimateLength(absEuclidianMatrix(s,d));
        absFiberMatrix(d,s) =  estimateLength(absEuclidianMatrix(d,s));
    end
end

% Adjacent matrix;
adjacentMatrix      = false(nNodes,nNodes);
for i = 1:(nNodes*(nNodes-1))
    [s,d] = ind2sub([nNodes,nNodes],i);
    if(distanceSpanMatrix(s,d)>=1e-4)
        adjacentMatrix(s,d)    = true;
        adjacentMatrix(d,s)    = true;
    end
end



% ==============================
% Export
% ==============================
save(exportMat...
    ,'absFiberMatrix' ...
    ,'absEuclidianMatrix' ...
    ,'distanceSpanMatrix' ...
    ,'latitudeDegreeArray' ...
    ,'longitudeDegreeArray' ...
    );

end
function fiber = estimateLength(euclidian)

if(euclidian<1e3*1e3)
    fiber = 1.5*euclidian;
elseif(euclidian<1.2e3*1e3)
    fiber = 1500*1e3;
else
    fiber = 1.25*euclidian;
end


end
function dist = fun_harversine(p1,p2)
%
% Description:
%   *Haversine formula* : calculating great-circle distance between two points over the earth's surface.
%  longitude: a geographic coordinate that specifies the east-west position of a point on the surface of the Earth
%  latitude :  ... south-north ...
%
% Example (calculating distance from Paris (N 48.8, E 2.4) to Shanghai (N 31.1,E 121.0) )
%     lat1 = 48.8*pi/180;
%     lat2 = 31.1*pi/180;
%     lon1 = 2.4*pi/180;
%     lon2 = 121*pi/180;
%     p1 = struct('latitude',lat1, ...
%         'longitude', lon1);
%     p2 = struct('latitude',lat2, ...
%         'longitude', lon2);
%     dist = fun_harversine(p1,p2);
%
% Reference: www.movable-type.co.uk/scripts/latlong.html
%
% ==============================
% Created by cao.chen;
% Date: 9th Jul., 2023
% ==============================

if(~(isa(p1,'struct')...
        &&isa(p2,'struct')))
    error('not defined');
end
lat1 = p1.latitude;
lat2 = p2.latitude;
lon1 = p1.longitude;
lon2 = p2.longitude;

dlat = abs(lat1 - lat2);
dlon = abs(lon1 - lon2);

Re   = 6367e3;                      % Earth radius; unit: km;
dist = 2*Re*asin(sqrt( sin(dlat/2)^2 + cos(lat1)*cos(lat2)* sin(dlon/2)^2 ));

end
