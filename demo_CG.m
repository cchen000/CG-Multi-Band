%% Demo: Column Generation (CG) for RWA in a 4-node network
%
%  This script solves the Routing and Wavelength Assignment (RWA) problem
%  using Column Generation (CG), an iterative decomposition method:
%    1. Restricted Master Problem (RMP): given a set of wavelength
%       configurations, decide how many times each configuration is used.
%    2. Pricing Sub-problem: find a new configuration with positive
%       reduced cost; if none exists, the current solution is optimal.
%  After CG converges with LP-relaxed RMP, re-solve RMP with integer
%  variables to obtain the final integer solution.
%
% Author: Cao Chen
% Date: 2026, Mar., 5th.
%
%% 1. Network parameters
% --- Resource limits ---
maxWavelengths  = 8;       % Total wavelength channels per fiber link (W)
maxTransceivers = 1000000; % Transceiver budget (set very large = unconstrained)

% --- Traffic demand matrix (4 nodes) ---
%  Entry (s,d) = proportion of total traffic from node s to node d.
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
Capacity124   = 100;
Capacity134   = 100;
Capacity14    = 100;

Capacity243   = 100;
Capacity213   = 100;
Capacity2143  = 50;

Capacity24    = 250;
Capacity2134  = 100;
Capacity214   = 50;

% --- Path name -> row index mapping ---
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

% String array of all path names (for display purposes)
pathStr = ["p124", 
    "p134", 
    "p14", 
    "p243", 
    "p213", 
    "p2143", 
    "p24", 
    "p2134", 
    "p214"];

% Solver options (suppress console output)
lpopts = optimoptions("linprog",    Display="off");
ipopts = optimoptions("intlinprog", Display="off");

%% 3. Initialize wavelength configuration pool
%  A "configuration" defines which paths are active on ONE wavelength.
%  Each configuration is described by:
%    T  (3x1) : carried capacity for each demand [d14, d23, d24]
%    isTransActive (9x1) : binary indicator — which of the 9 paths are ON
%    noTrans : number of active transceivers in this configuration
%
%  We seed the pool with 3 trivial single-path configurations.

% Config 1: only path p124 active -> serves demand (1,4) with 100 Gbps
T1             = [100;0;0];
isTransActive1 = [1;0;0;0;0;0;0;0;0];
noTrans1       = sum(isTransActive1);

% Config 2: only path p243 active -> serves demand (2,3) with 100 Gbps
T2             = [0;100;0];
isTransActive2 = [0;0;0;1;0;0;0;0;0];
noTrans2       = sum(isTransActive2);

% Config 3: only path p24 active  -> serves demand (2,4) with 250 Gbps
T3             = [0;0;250]; 
isTransActive3 = [0;0;0;0;0;0;1;0;0];
noTrans3       = sum(isTransActive3);

% Aggregate initial columns into matrices (each column = one config)
tranCapacityArray = [T1,T2,T3];         % 3 x nConfigs
isTransActiveArray= [isTransActive1, isTransActive2, isTransActive3]; % 9 x nConfigs

noConfigs         = size(tranCapacityArray,1);

NumTransArray     = [noTrans1, noTrans2, noTrans3]; % 1 x nConfigs

%% 4. Column Generation loop
%  Iterate between RMP and Pricing until no improving column exists.
reducedCost       = inf;             % initialize to enter the loop
tranCapacityArray = [T1,T2,T3];      % reset column pool
sols              = struct('T',[], 'a', [], 'x',[]);

while reducedCost >= 1e-6
    % -------------------------------------------------------
    % 4a. Append the new configuration found by Pricing
    %     (first iteration: sols fields are empty, so nothing is added)
    % -------------------------------------------------------
    tranCapacityArray  = [tranCapacityArray,  sols.T' ];
    NumTransArray      = [NumTransArray,      sols.a'];
    isTransActiveArray = [isTransActiveArray, sols.x' ];
    
    % -------------------------------------------------------
    % 4b. Solve Restricted Master Problem (LP-relaxed)
    % -------------------------------------------------------
    %  Decision variable z(k): how many wavelengths use configuration k
    %  (continuous relaxation -> fractional z allowed)
    rwaprob = ...
          optimproblem('Description','Max Throughput - RMP',...
          'ObjectiveSense', 'maximize');
    z   = optimvar('z',size(tranCapacityArray,2), ...
        'LowerBound', 0, ...
        'Type', 'continuous');
    TN  = optimvar('TN', 1, 'LowerBound', 0, 'Type', 'continuous');
    
    % Demand-capacity coupling: carried capacity >= proportional share
    rwaprob.Constraints.DemandCapcity14  = ...
        normalized_demand_matrix(1,4)*TN <= tranCapacityArray(1,:)*z;
    rwaprob.Constraints.DemandCapcity23  = ...
        normalized_demand_matrix(2,3)*TN <= tranCapacityArray(2,:)*z;
    rwaprob.Constraints.DemandCapcity24  = ...
        normalized_demand_matrix(2,4)*TN <= tranCapacityArray(3,:)*z;

    % Total wavelength usage <= W
    rwaprob.Constraints.Wavelength = sum(z)<=maxWavelengths;

    % Transceiver budget
    rwaprob.Constraints.TransceiverCost = NumTransArray*z<=maxTransceivers;

    % Objective: maximize network throughput
    rwaprob.Objective.netwThroughput = TN;

    [valueRMP, maxThroughput, exitflag,~,lambda_RMP] = solve(rwaprob,'options',lpopts);

    % --- Extract dual variables (shadow prices) ---
    %  MATLAB returns Lagrange multipliers with opposite sign to the
    %  standard dual variable convention, so we flip them here.
    signOfMatlabSolver = -1;
    lambda_RMP.Constraints.DemandCapcity14 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity14;
    lambda_RMP.Constraints.DemandCapcity23 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity23;
    lambda_RMP.Constraints.DemandCapcity24 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity24;
    lambda_RMP.Constraints.Wavelength      = signOfMatlabSolver * lambda_RMP.Constraints.Wavelength;
    lambda_RMP.Constraints.TransceiverCost = signOfMatlabSolver * lambda_RMP.Constraints.TransceiverCost;

    % Print RMP solution and dual values
    fprintf('(Relaxed RMP) Maximal network throughput %g [Gbps]\n',valueRMP.TN);
    fprintf("\t\tz=[%s]\n", strjoin(string(valueRMP.z'), ", "));
    vals = structfun(@(x) x, lambda_RMP.Constraints);
    names = fieldnames(lambda_RMP.Constraints);
    fprintf("\t\t[%s]", strjoin(compose("%s=%.4f", string(names), vals), ", "));
    fprintf('\n');

    % -------------------------------------------------------
    % 4c. Solve Pricing Sub-problem
    % -------------------------------------------------------
    %  Find a NEW wavelength configuration (x, T, a) that has the
    %  largest reduced cost.  If reduced cost <= 0, no improving
    %  column exists and CG terminates.
    %
    %  Variables:
    %    x(p) in {0,1} : whether path p is active in this configuration
    %    T(d)          : capacity delivered to demand d
    %    a             : number of transceivers used
    subproblem = optimproblem('Description','New Configuration - Pricing',...
        'ObjectiveSense', 'maximize');
    T  = optimvar('T', ["14","23","24"], 'LowerBound', 0, 'Type', 'continuous');
    x  = optimvar('x', {'p124', 'p134','p14', ...
        'p243','p213','p2143', ...
        'p24','p2134','p214'}, ...
        'LowerBound', 0, 'Type', 'integer');
    a  = optimvar('a', 'LowerBound', 0, 'Type', 'integer');

    % Wavelength clash-free: at most one path per link in this config
    subproblem.Constraints.WavelengthNonOverlapping12  = ...  % Link 1->2
        x('p124') <= 1;
    subproblem.Constraints.WavelengthNonOverlapping13 = ...   % Link 1->3
        x('p134') + x('p213') + x('p2134')  <= 1;
    subproblem.Constraints.WavelengthNonOverlapping14 = ...   % Link 1->4
        x('p14') + x('p2143') + x('p214')  <= 1;
    subproblem.Constraints.WavelengthNonOverlapping24 = ...   % Link 2->4
        x('p124') + x('p243') + x('p24')   <= 1;
    subproblem.Constraints.WavelengthNonOverlapping34 = ...   % Link 3->4
        x('p134') + x('p2134')             <= 1;
    subproblem.Constraints.WavelengthNonOverlapping21 = ...   % Link 2->1
        x('p213') + x('p2143')  + x('p2134') + x('p214') <= 1;
    subproblem.Constraints.WavelengthNonOverlapping43 = ...   % Link 4->3
        x('p243') + x('p2143')             <= 1;

    % Capacity delivered per demand = sum of active path capacities
    subproblem.Constraints.DemandCapcity14  = ...
        T('14') == sum(Capacity124 * x('p124') + Capacity134 * x('p134') + Capacity14   * x('p14'));
    subproblem.Constraints.DemandCapcity23  = ...
        T('23') == sum(Capacity243 * x('p243') + Capacity213 * x('p213') + Capacity2143 * x('p2143'));
    subproblem.Constraints.DemandCapcity24  = ...
        T('24') == sum(Capacity24  * x('p24')  + Capacity2134* x('p2134')+ Capacity214  * x('p214'));

    % Transceiver count = total active paths
    subproblem.Constraints.CostLimit = ...
        sum(x)==a;

    % Reduced cost = dual(demand) * T - dual(transceiver) * a - dual(wavelength)
    subproblem.Objective.ReducedCost =  ...
        + lambda_RMP.Constraints.DemandCapcity14*T('14') ...
        + lambda_RMP.Constraints.DemandCapcity23*T('23') ...
        + lambda_RMP.Constraints.DemandCapcity24*T('24') ...
        - lambda_RMP.Constraints.TransceiverCost*a ...
        - lambda_RMP.Constraints.Wavelength ...
        ;

    [sols,reducedCost,exitflag,~,~] = solve(subproblem,'options',ipopts);

    % --- Check termination ---
    if reducedCost <= 1e-6
        fprintf('(Pricing) Reduced cost = %g',reducedCost);
        fprintf('\t\tNo improving column -> CG converged\n\n');
        break;
    else
        fprintf('(Pricing) Reduced cost = %g',reducedCost);
        fprintf('\t\tAdd configuration {%s}\n\n',strjoin(pathStr(logical(sols.x)),','));
    end
end


%% 5. Re-solve RMP with integer variables
%  The CG loop used LP relaxation; now fix z to integers for a valid
%  wavelength assignment.

z.Type = 'integer';
[valueRMP, maxThroughput, exitflag,~,lambda_RMP] = solve(rwaprob,'options',ipopts);

fprintf('Optimal solutions:\n');
fprintf(' (Integer RMP) Maximal network throughput %g [Gbps]\n',valueRMP.TN);
for i = 1:length(valueRMP.z)
    nConfigs = round(valueRMP.z(i));
    if nConfigs >= 1
        % Find which paths are active in configuration i
        activeIdx = cellfun(@(p) isTransActiveArray(pathIdx(p), i) > 0, pathStr);
        pathNames   = strjoin(string(pathStr(activeIdx)), ", ");
        fprintf(' configurations {%s} in %d times \n', pathNames, nConfigs);
    end
end

%% 6. Map configurations to actual wavelengths and visualize
%  Each configuration k occupies z(k) consecutive wavelength slots.
%  We walk through configurations and fill the (link x wavelength) matrix.

wStart = 1;  % next available wavelength slot
link12 = zeros(1,maxWavelengths);
link13 = zeros(1,maxWavelengths);
link14 = zeros(1,maxWavelengths);
link24 = zeros(1,maxWavelengths);
link34 = zeros(1,maxWavelengths);
link21 = zeros(1,maxWavelengths);
link43 = zeros(1,maxWavelengths);

for i = 1:length(valueRMP.z)
    nConfigs = round(valueRMP.z(i));
    if nConfigs >= 1
        wRange = wStart:(wStart+nConfigs-1);  % wavelength slots for config i

        % For each link, mark occupied if any traversing path is active
        link12(wRange) = isTransActiveArray(pathIdx('p124'),i);

        link13(wRange) = isTransActiveArray(pathIdx('p134'),i)  +  isTransActiveArray(pathIdx('p213'),i)  + isTransActiveArray(pathIdx('p2134'),i);

        link14(wRange) = isTransActiveArray(pathIdx('p14'),i)  +  isTransActiveArray(pathIdx('p2143'),i)  + isTransActiveArray(pathIdx('p214'),i);

        link24(wRange) = isTransActiveArray(pathIdx('p124'),i)  +  isTransActiveArray(pathIdx('p243'),i)  + isTransActiveArray(pathIdx('p24'),i) ;

        link34(wRange) = isTransActiveArray(pathIdx('p134'),i)  +  isTransActiveArray(pathIdx('p2134'),i);

        link21(wRange) = isTransActiveArray(pathIdx('p213'),i)  ...
            +  isTransActiveArray(pathIdx('p2143'),i) +   isTransActiveArray(pathIdx('p2134'),i) +   isTransActiveArray(pathIdx('p214'),i);

        link43(wRange) = isTransActiveArray(pathIdx('p243'),i)  +  isTransActiveArray(pathIdx('p2143'),i) ;
    end
    wStart = wStart + nConfigs;  % advance to next free slot
end

% --- Assemble into a (links x wavelengths) matrix and plot ---
networkSpectrum = [link12;link13;link14;link24;link34;link21;link43];

%% 7. Plot wavelength usage heatmap
imagesc(networkSpectrum);  % yellow = occupied, blue = free
colorbar;
xlabel("wavelength channel");
ylabel("link ID");
title('Link-Wavelength Occupation Map');
xlim([1-0.5 8+0.5]);
ylim([1-0.5 7+0.5]);
