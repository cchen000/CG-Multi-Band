%% Demo: ILP-based Routing and Wavelength Assignment (RWA) in a 4-node network
%
%  This script solves the RWA problem using Integer Linear Programming (ILP).
%  Goal: Maximize total network throughput by jointly deciding:
%        (1) which candidate path to use for each demand,
%        (2) which wavelength to assign on each path,
%  subject to wavelength non-overlapping (clash-free) constraints on every link.
%
% Author: Cao Chen
% Date: 2026, Mar., 5th.
%
%% 1. Network parameters
% 
% --- Resource limits ---
maxWavelengths  = 8;       % Total number of wavelength channels per fiber link (W)
maxTransceivers = 1000000; % Total transceiver budget (set very large = unconstrained)

W               = maxWavelengths;
A               = maxTransceivers;

% --- Traffic demand matrix (4 nodes) ---
%  Entry (s,d) gives the *proportion* of total traffic from node s to node d.
%  Three non-zero demands: (1->4), (2->3), (2->4)
normalized_demand_matrix = [
  0 0 0 1/3
  0 0 1/3 1/3
  0 0 0 0
  0 0 0 0  
];

% ==============================
% 2. Candidate paths (3 paths per demand, 9 paths total)
% ==============================
%  Topology (4 nodes, 7 directed links):
%       1 ---8--- 2
%       |\        |
%       4  25     2
%       |    \    |
%       3 ---9--- 4
%  (numbers on edges = link distances)

% --- Demand 1->4: three candidate paths ---
Path124     = [1,2,4];    % 1 -> 2 -> 4
Path134     = [1,3,4];    % 1 -> 3 -> 4
Path14      = [1,4];      % 1 -> 4 (direct)

% --- Demand 2->3: three candidate paths ---
Path243     = [2,4,3];    % 2 -> 4 -> 3
Path213     = [2,1,3];    % 2 -> 1 -> 3
Path2143    = [2,1,4,3];  % 2 -> 1 -> 4 -> 3

% --- Demand 2->4: three candidate paths ---
Path24      = [2,4];      % 2 -> 4 (direct)
Path2134    = [2,1,3,4];  % 2 -> 1 -> 3 -> 4
Path214     = [2,1,4];    % 2 -> 1 -> 4

% --- Route length (sum of link distances along each path) ---
RouteLength124  = 10;
RouteLength134  = 13;
RouteLength14   = 25;
RouteLength243  = 11;
RouteLength213  = 12;
RouteLength2143 = 42;
RouteLength24   = 2;
RouteLength2134 = 21;
RouteLength214  = 33;

% --- Capacity per lightpath on each route [Gbps] ---
%  (shorter / better-quality routes may carry higher capacity)
Capacity124   = 100;
Capacity134   = 100;
Capacity14    = 100;

Capacity243   = 100;
Capacity213   = 100;
Capacity2143  = 50;

Capacity24    = 250;
Capacity2134  = 100;
Capacity214   = 50;

% --- Path name -> row index mapping (for reading solution later) ---
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


%% 3. Build ILP model
rwaprob = optimproblem('Description','Max Throughput - ILP');

% --- Decision variable x(path, wavelength) ---
%  x(p, w) = 1  if path p is assigned wavelength w
%  x(p, w) = 0  otherwise
%  Dimensions: 9 paths x W wavelengths
x = optimvar('x', {'p124', 'p134','p14', ...
    'p243','p213','p2143', ...
    'p24','p2134','p214'}, ...
    maxWavelengths, ...
    'LowerBound', 0, 'Type', 'integer');

% --- Auxiliary variables ---
%  T(d)  = total carried capacity for each demand d
%  TN    = total network throughput (to be maximized)
T  = optimvar('T', ["14","23","24"], 'LowerBound', 0, 'Type', 'continuous');
TN = optimvar('TN', 1, 'LowerBound', 0, 'Type', 'continuous');


% ==============================
% Constraint 1: Wavelength clash-free on every link
% ==============================
%  On each directed link, at most ONE path may use the same wavelength.
%  For each link e and wavelength w:  sum of x(p,w) over all paths p
%  that traverse link e  <=  1.

rwaprob.Constraints.WavelengthNonOverlapping12  = ...  % Link 1->2: only p124 uses it
    x('p124',1:W) <= 1;

rwaprob.Constraints.WavelengthNonOverlapping13 = ...   % Link 1->3: p134, p213, p2134
    x('p134',1:W) + x('p213',1:W) + x('p2134',1:W)  <= 1;

rwaprob.Constraints.WavelengthNonOverlapping14 = ...   % Link 1->4: p14, p2143, p214
    x('p14',1:W) + x('p2143',1:W) + x('p214',1:W)  <= 1;

rwaprob.Constraints.WavelengthNonOverlapping24 = ...   % Link 2->4: p124, p243, p24
    x('p124',1:W) + x('p243',1:W) + x('p24',1:W)  <= 1;

rwaprob.Constraints.WavelengthNonOverlapping34 = ...   % Link 3->4: p134, p2134
    x('p134',1:W) + x('p2134',1:W)                <= 1;

rwaprob.Constraints.WavelengthNonOverlapping21 = ...   % Link 2->1: p213, p2143, p2134, p214
    x('p213',1:W) + x('p2143',1:W)  + x('p2134',1:W) + x('p214',1:W)  <= 1;

rwaprob.Constraints.WavelengthNonOverlapping43 = ...   % Link 4->3: p243, p2143
    x('p243',1:W) + x('p2143',1:W)               <= 1;


% ==============================
% Constraint 2: Transceiver budget
% ==============================
%  Total activated lightpaths must not exceed available transceivers.
rwaprob.Constraints.CostLimit = ...
    sum(sum(x))<=A;

% ==============================
% Constraint 3: Demand-capacity linkage
% ==============================
%  T(d) = total capacity carried for demand d
%       = sum over all (path, wavelength) assigned to demand d of Capacity(path).
rwaprob.Constraints.DemandCapcity14  = ...
    T('14') == sum(sum(Capacity124 * x('p124',:) + Capacity134 * x('p134',:) + Capacity14 * x('p14',:)));
rwaprob.Constraints.DemandCapcity23  = ...
    T('23') == sum(sum(Capacity243 * x('p243',:) + Capacity213 * x('p213',:) + Capacity2143 * x('p2143',:)));
rwaprob.Constraints.DemandCapcity24  = ...
    T('24') == sum(sum(Capacity24  * x('p24',:)  + Capacity2134* x('p2134',:) + Capacity214 * x('p214',:)));

% ==============================
% Constraint 4: Proportional fairness (bottleneck formulation)
% ==============================
%  Network throughput TN is limited by the "worst" demand ratio:
%  TN <= T(d) / normalized_demand(d)  for each demand d.
%  Maximizing TN ensures all demands scale proportionally.
rwaprob.Constraints.DemandNetwork14 = ...
    TN <= 1/normalized_demand_matrix(1,4) * T('14');
rwaprob.Constraints.DemandNetwork23 = ...
    TN <= 1/normalized_demand_matrix(2,3) * T('23');
rwaprob.Constraints.DemandNetwork24 = ...
    TN <= 1/normalized_demand_matrix(2,4) * T('24');

% --- Objective: maximize TN (negate because MATLAB minimizes by default) ---
rwaprob.Objective.netwThroughput = -TN;

%% 4. Solve
lpopts = optimoptions("linprog",Display="off");
ipopts = optimoptions("intlinprog",Display="off");

[values,nLogs,exitflag,~,lambda] = solve(rwaprob,'options',ipopts);
fprintf('Maximal network throughput %g [Gbps]\n',values.TN);


%% 5. Visualize wavelength usage on each link
%  For every directed link, sum up x(p,w) over all paths that traverse it.
%  Result: 1 = wavelength occupied, 0 = free.

% Link 1->2
link12 =  values.x(pathIdx('p124'),:);

% Link 1->3
link13 =  values.x(pathIdx('p134'),:) + values.x(pathIdx('p213'),:) + values.x(pathIdx('p2134'),:);

% Link 1->4
link14 =  values.x(pathIdx('p14'),:) + values.x(pathIdx('p2143'),:) + values.x(pathIdx('p214'),:);

% Link 2->4
link24 =  values.x(pathIdx('p124'),:) + values.x(pathIdx('p243'),:) + values.x(pathIdx('p24'),:);

% Link 3->4
link34 =  values.x(pathIdx('p134'),:) + values.x(pathIdx('p2134'),:);

% Link 2->1
link21 =  values.x(pathIdx('p213'),:) + values.x(pathIdx('p2143'),:)  + values.x(pathIdx('p2134'),:) ...
    + values.x(pathIdx('p214'),:);

% Link 4->3
link43 =  values.x(pathIdx('p243'),:) + values.x(pathIdx('p2143'),:);

% --- Assemble into a (links x wavelengths) matrix and plot ---
networkSpectrum = [link12;link13;link14;link24;link34;link21;link43];

imagesc(networkSpectrum);  % heatmap: yellow=occupied, blue=free
colorbar;
xlabel("wavelength channel");
ylabel("link ID");
xlim([1-0.5 8+0.5]);
ylim([1-0.5 7+0.5]);