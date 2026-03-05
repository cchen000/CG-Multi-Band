%% Demo implementation of CG in 4-node network
% 
% Author: Cao Chen
% Date: 2026, Mar., 5th.
% 
%% Loading basic parameters
maxWavelengths  = 8; % Maximum available wavelengths, W
maxTransceivers = 1000000; % Maximum availble transceivers, A

W               = maxWavelengths;
A               = maxTransceivers;

normalized_demand_matrix = [
  0 0 0 1/3
  0 0 1/3 1/3
  0 0 0 0
  0 0 0 0  
]; % Network traffic demand distribution;

% ==============================
% 3 candidate paths for 3 demands;
% ==============================

% Optical Network data;
% OpticalNetwork = zeros(4,4);
% OpticalNetwork(1,2) = 8;
% OpticalNetwork(2,1) = 8;
% OpticalNetwork(1,3) = 4;
% OpticalNetwork(1,4) = 25;
% OpticalNetwork(2,4) = 2;
% OpticalNetwork(3,4) = 9;
% OpticalNetwork(4,3) = 9;

% Path
Path124     = [1,2,4];
Path134     = [1,3,4];
Path14      = [1,4];

Path243     = [2,4,3];
Path213     = [2,1,3];
Path2143    = [2,1,4,3];

Path24      = [2,4];
Path2134    = [2,1,3,4];
Path214     = [2,1,4];

% RouteLength;
RouteLength124  = 10;
RouteLength134  = 13;
RouteLength14   = 25;
RouteLength243  = 11;
RouteLength213  = 12;
RouteLength2143 = 42;
RouteLength24   = 2;
RouteLength2134 = 21;
RouteLength214  = 33;

% Capacity;
Capacity124   = 100;
Capacity134   = 100;
Capacity14    = 100;

Capacity243   = 100;
Capacity213   = 100;
Capacity2143  = 50;

Capacity24    = 250;
Capacity2134  = 100;
Capacity214   = 50;

% PathIdx
pathIdx = containers.Map('KeyType', 'char', 'ValueType', 'int32');
pathIdx('p124') = 1;
pathIdx('p134') = 2;
pathIdx('p14')  = 3;
pathIdx('p243') = 4;
pathIdx('p213') = 5;
pathIdx('p2143')= 6;
pathIdx('p24')  = 7;
pathIdx('p2134')= 8;
pathIdx('p214') = 9;


%% Model
rwaprob = optimproblem('Description','Max Throughput - ILP');
% Create variables representing the number of each pattern used
x = optimvar('x', {'p124', 'p134','p14', ...
    'p243','p213','p2143', ...
    'p24','p2134','p214'}, ...
    maxWavelengths, ...
    'LowerBound', 0, 'Type', 'integer');


T  = optimvar('T', ["14","23","24"], 'LowerBound', 0, 'Type', 'continuous');
TN = optimvar('TN', 1, 'LowerBound', 0, 'Type', 'continuous');


% ==============================
% Wavelength non-overlapping for links;
% ==============================

rwaprob.Constraints.WavelengthNonOverlapping12  = ...
    x('p124',1:W) <= 1; % Non-overlapping for LINK12

rwaprob.Constraints.WavelengthNonOverlapping13 = ...
    x('p134',1:W) + x('p213',1:W) + x('p2134',1:W)  <= 1;%  LINK13

rwaprob.Constraints.WavelengthNonOverlapping14 = ...
    x('p14',1:W) + x('p2143',1:W) + x('p214',1:W)  <= 1;% LINK14

rwaprob.Constraints.WavelengthNonOverlapping24 = ...
    x('p124',1:W) + x('p243',1:W) + x('p24',1:W)  <= 1; % LINK24

rwaprob.Constraints.WavelengthNonOverlapping34 = ...
    x('p134',1:W) + x('p2134',1:W)                <= 1; %  LINK34

rwaprob.Constraints.WavelengthNonOverlapping21 = ...% LINK21
    x('p213',1:W) + x('p2143',1:W)  + x('p2134',1:W) + x('p214',1:W)  <= 1;

rwaprob.Constraints.WavelengthNonOverlapping43 = ...
    x('p243',1:W) + x('p2143',1:W)               <= 1;% LINK43


% ==============================
% max. Transsceivers;
% ==============================
rwaprob.Constraints.CostLimit = ...
    sum(sum(x))<=A;


% ==============================
% Demand-capacity constraints
% ==============================
rwaprob.Constraints.DemandCapcity14  = ...
    T('14') == sum(sum(Capacity124 * x('p124',:) + Capacity134 * x('p134',:) + Capacity14 * x('p14',:)));
rwaprob.Constraints.DemandCapcity23  = ...
    T('23') == sum(sum(Capacity243 * x('p243',:) + Capacity213 * x('p213',:) + Capacity2143 * x('p2143',:)));
rwaprob.Constraints.DemandCapcity24  = ...
    T('24') == sum(sum(Capacity24  * x('p24',:)  + Capacity2134* x('p2134',:) + Capacity214 * x('p214',:)));

rwaprob.Constraints.DemandNetwork14 = ...
    TN <= 1/normalized_demand_matrix(1,4) * T('14');
rwaprob.Constraints.DemandNetwork23 = ...
    TN <= 1/normalized_demand_matrix(2,3) * T('23');
rwaprob.Constraints.DemandNetwork24 = ...
    TN <= 1/normalized_demand_matrix(2,4) * T('24');


% ==============================
% Transciever cost;
% ==============================
rwaprob.Constraints.TransceiverCost = sum(sum(sum(x)))<=maxTransceivers;

% The objective is the number of logs used
rwaprob.Objective.netwThroughput = -TN;

lpopts = optimoptions("linprog",Display="off");
ipopts = optimoptions("intlinprog",Display="off");

[values,nLogs,exitflag,~,lambda] = solve(rwaprob,'options',ipopts);
fprintf('Maximal network throughput %g [Gbps]\n',values.TN);


%% Plot

c = 3;
% LINK12
link12 =  values.x(pathIdx('p124'),:); % 'p124'

% LINK13
link13 =  values.x(pathIdx('p134'),:) + values.x(pathIdx('p213'),:) + values.x(pathIdx('p2134'),:);

% LINK14
link14 =  values.x(pathIdx('p14'),:) + values.x(pathIdx('p2143'),:) + values.x(pathIdx('p214'),:);

% LINK24
link24 =  values.x(pathIdx('p124'),:) + values.x(pathIdx('p243'),:) + values.x(pathIdx('p24'),:);

% LINK34
link34 =  values.x(pathIdx('p134'),:) + values.x(pathIdx('p2134'),:);

% LINK21
link21 =  values.x(pathIdx('p213'),:) + values.x(pathIdx('p2143'),:)  + values.x(pathIdx('p2134'),:) ...
    + values.x(pathIdx('p214'),:);

% LINK43
link43 =  values.x(pathIdx('p243'),:) + values.x(pathIdx('p2143'),:);

networkSpectrum = [link12;link13;link14;link24;link34;link21;link43];


% % 
% 
imagesc(networkSpectrum);
colorbar;
xlabel("wavelength channel");
ylabel("link ID");
xlim([1-0.5 8+0.5]);
ylim([1-0.5 7+0.5]);