%% Demo implementation of CG in 4-node network
% 
% Author: Cao Chen
% Date: 2026, Mar., 5th.
% 
%% Loading parameters
maxWavelengths  = 8; % Maximum available wavelengths, W in the paper;
maxTransceivers = 1000000; % Maximum availble transceivers, A  in the paper;

normalized_demand_matrix = [
  0 0 0 1/3
  0 0 1/3 1/3
  0 0 0 0
  0 0 0 0  
]; % Network traffic demand distribution;

% ==============================
% 3 candidate paths for 3 demands;
% ==============================

% Optical Network Data
% OpticalNetwork = zeros(4,4);
% OpticalNetwork(1,2) = 8;
% OpticalNetwork(2,1) = 8;
% OpticalNetwork(1,3) = 4;
% OpticalNetwork(1,4) = 25;
% OpticalNetwork(2,4) = 2;
% OpticalNetwork(3,4) = 9;
% OpticalNetwork(4,3) = 9;

% Path for demands (1,4), (2,3), (2,4)
Path124     = [1,2,4];
Path134     = [1,3,4];
Path14      = [1,4];

Path243     = [2,4,3];
Path213     = [2,1,3];
Path2143    = [2,1,4,3];

Path24      = [2,4];
Path2134    = [2,1,3,4];
Path214     = [2,1,4];

% Path Span Number
RouteLength124  = 10;
RouteLength134  = 13;
RouteLength14   = 25;
RouteLength243  = 11;
RouteLength213  = 12;
RouteLength2143 = 42;
RouteLength24   = 2;
RouteLength2134 = 21;
RouteLength214  = 33;

% Path Capacity
Capacity124   = 100;
Capacity134   = 100;
Capacity14    = 100;

Capacity243   = 100;
Capacity213   = 100;
Capacity2143  = 50;

Capacity24    = 250;
Capacity2134  = 100;
Capacity214   = 50;

% Path Index;
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

pathStr = ["p124", 
    "p134", 
    "p14", 
    "p243", 
    "p213", 
    "p2143", 
    "p24", 
    "p2134", 
    "p214"];

lpopts = optimoptions("linprog",    Display="off");
ipopts = optimoptions("intlinprog", Display="off");

%% Initialize wavelength configuration set

% Config1 : 'p124'
T1             = [100;0;0]; % demands (1,4), (2,3), (2,4)
isTransActive1 = [1;0;0;0;0;0;0;0;0];
noTrans1       = sum(isTransActive1);

% Config2  : 'p243'
T2             = [0;100;0];
isTransActive2 = [0;0;0;1;0;0;0;0;0];
noTrans2       = sum(isTransActive2);

% Config3  : 'p24'
T3             = [0;0;250]; 
isTransActive3 = [0;0;0;0;0;0;1;0;0];
noTrans3       = sum(isTransActive3);

tranCapacityArray = [T1,T2,T3];
isTransActiveArray= [isTransActive1, isTransActive2, isTransActive3];

noConfigs         = size(tranCapacityArray,1);

NumTransArray     = [noTrans1, noTrans2, noTrans3];

%% while - loop
reducedCost       = inf;
tranCapacityArray = [T1,T2,T3];
sols              = struct('T',[], 'a', [], 'x',[]);
while reducedCost>=1e-6
    %% Add new column into the set
    tranCapacityArray  = [tranCapacityArray,  sols.T' ];
    NumTransArray      = [NumTransArray,      sols.a'];
    isTransActiveArray = [isTransActiveArray, sols.x' ];
    
    %% solve Restricted Master Problem with relaxing variables
    % Create variables representing the number of each pattern used
    rwaprob = ...
          optimproblem('Description','Max Throughput - RMP',...
          'ObjectiveSense', 'maximize');
    z   = optimvar('z',size(tranCapacityArray,2), ...
        'LowerBound', 0, ...
        'Type', 'continuous');
    TN  = optimvar('TN', 1, 'LowerBound', 0, 'Type', 'continuous');
    
    rwaprob.Constraints.DemandCapcity14  = ...
        normalized_demand_matrix(1,4)*TN <= tranCapacityArray(1,:)*z;
    rwaprob.Constraints.DemandCapcity23  = ...
        normalized_demand_matrix(2,3)*TN <= tranCapacityArray(2,:)*z;
    rwaprob.Constraints.DemandCapcity24  = ...
        normalized_demand_matrix(2,4)*TN <= tranCapacityArray(3,:)*z;
    rwaprob.Constraints.Wavelength = sum(z)<=maxWavelengths;
    rwaprob.Constraints.TransceiverCost = NumTransArray*z<=maxTransceivers;

    % The objective is the number of logs used
    rwaprob.Objective.netwThroughput = TN;

    [valueRMP, maxThroughput, exitflag,~,lambda_RMP] = solve(rwaprob,'options',lpopts);
    % % NOTE: FOR MATLAB, LAMBDA is LAGRANGE MULTIPLIXER;
    % % We only need the dual variable, which is the opposite of the lagrange multiplexer;
    % % So, we should convert it by ourselves
    % % This will be used in dual objective value;
    signOfMatlabSolver = -1; % to all lambda;
    lambda_RMP.Constraints.DemandCapcity14 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity14;
    lambda_RMP.Constraints.DemandCapcity23 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity23;
    lambda_RMP.Constraints.DemandCapcity24 = signOfMatlabSolver * lambda_RMP.Constraints.DemandCapcity24;
    lambda_RMP.Constraints.Wavelength      = signOfMatlabSolver * lambda_RMP.Constraints.Wavelength;
    lambda_RMP.Constraints.TransceiverCost = signOfMatlabSolver * lambda_RMP.Constraints.TransceiverCost;

    fprintf('(Relaxed RMP) Maximal network throughput %g [Gbps]\n',valueRMP.TN);
    fprintf("\t\tz=[%s]\n", strjoin(string(valueRMP.z'), ", "));
    vals = structfun(@(x) x, lambda_RMP.Constraints);
    names = fieldnames(lambda_RMP.Constraints);
    fprintf("\t\t[%s]", strjoin(compose("%s=%.4f", string(names), vals), ", "));
    fprintf('\n');

    %% solve Pricing problem and Choose a column
    subproblem = optimproblem('Description','New Configuration - Pricing',...
        'ObjectiveSense', 'maximize');
    T  = optimvar('T', ["14","23","24"], 'LowerBound', 0, 'Type', 'continuous');
    x  = optimvar('x', {'p124', 'p134','p14', ...
        'p243','p213','p2143', ...
        'p24','p2134','p214'}, ...
        'LowerBound', 0, 'Type', 'integer');
    a  = optimvar('a', 'LowerBound', 0, 'Type', 'integer');

    subproblem.Constraints.WavelengthNonOverlapping12  = ...
        x('p124') <= 1; % Non-overlapping for LINK12
    subproblem.Constraints.WavelengthNonOverlapping13 = ...
        x('p134') + x('p213') + x('p2134')  <= 1;%  LINK13
    subproblem.Constraints.WavelengthNonOverlapping14 = ...
        x('p14') + x('p2143') + x('p214')  <= 1;% LINK14
    subproblem.Constraints.WavelengthNonOverlapping24 = ...
        x('p124') + x('p243') + x('p24')   <= 1; % LINK24
    subproblem.Constraints.WavelengthNonOverlapping34 = ...
        x('p134') + x('p2134')             <= 1; %  LINK34
    subproblem.Constraints.WavelengthNonOverlapping21 = ...% LINK21
        x('p213') + x('p2143')  + x('p2134') + x('p214') <= 1;
    subproblem.Constraints.WavelengthNonOverlapping43 = ...
        x('p243') + x('p2143')             <= 1;% LINK43

    subproblem.Constraints.DemandCapcity14  = ...
        T('14') == sum(Capacity124 * x('p124') + Capacity134 * x('p134') + Capacity14   * x('p14'));
    subproblem.Constraints.DemandCapcity23  = ...
        T('23') == sum(Capacity243 * x('p243') + Capacity213 * x('p213') + Capacity2143 * x('p2143'));
    subproblem.Constraints.DemandCapcity24  = ...
        T('24') == sum(Capacity24  * x('p24')  + Capacity2134* x('p2134')+ Capacity214  * x('p214'));

    subproblem.Constraints.CostLimit = ...
        sum(x)==a;

    subproblem.Objective.ReducedCost =  ...
        + lambda_RMP.Constraints.DemandCapcity14*T('14') ...
        + lambda_RMP.Constraints.DemandCapcity23*T('23') ...
        + lambda_RMP.Constraints.DemandCapcity24*T('24') ...
        - lambda_RMP.Constraints.TransceiverCost*a ...
        - lambda_RMP.Constraints.Wavelength ...
        ;
    [sols,reducedCost,exitflag,~,~] = solve(subproblem,'options',ipopts);
    if reducedCost <= 1e-6
        fprintf('(Prcing)Reduced cost =%g',reducedCost);
        fprintf('\t\tskip\n\n');
        break;
    else
        fprintf('(Prcing)Reduced cost =%g',reducedCost);
        fprintf('\t\tAdd configuration {%s}\n\n',strjoin(pathStr(logical(sols.x)),','));
    end
end


%% Solve RMP with integer variables

z.Type = 'integer';
[valueRMP, maxThroughput, exitflag,~,lambda_RMP] = solve(rwaprob,'options',ipopts);

fprintf('Optimal solutions:\n');
fprintf(' (Integer RMP) Maximal network throughput %g [Gbps]\n',valueRMP.TN);
for i = 1:length(valueRMP.z)
    nConfigs = round(valueRMP.z(i));
    if nConfigs >= 1
        activeIdx = cellfun(@(p) isTransActiveArray(pathIdx(p), i) > 0, pathStr);
        pathNames   = strjoin(string(pathStr(activeIdx)), ", ");
        fprintf(' configurations {%s} in %d times \n', pathNames, nConfigs);
    end
end

%% Assign a wavelength to each configuration
wStart = 1;
i      = 1;
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
        % LINK12
        link12((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p124'),i);

        % LINK13
        link13((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p134'),i)  +  isTransActiveArray(pathIdx('p213'),i)  + isTransActiveArray(pathIdx('p2134'),i);

        % LINK14
        link14((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p14'),i)  +  isTransActiveArray(pathIdx('p2143'),i)  + isTransActiveArray(pathIdx('p214'),i);

        % LINK24
        link24((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p124'),i)  +  isTransActiveArray(pathIdx('p243'),i)  + isTransActiveArray(pathIdx('p24'),i) ;

        % LINK34
        link34((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p134'),i)  +  isTransActiveArray(pathIdx('p2134'),i);

        % LINK21
        link21((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p213'),i)  ...
            +  isTransActiveArray(pathIdx('p2143'),i) +   isTransActiveArray(pathIdx('p2134'),i) +   isTransActiveArray(pathIdx('p214'),i);

        % LINK43
        link43((wStart: (wStart+nConfigs-1))) = isTransActiveArray(pathIdx('p243'),i)  +  isTransActiveArray(pathIdx('p2143'),i) ;

    end
    wStart = wStart + nConfigs; % update new beginning
end

networkSpectrum = [link12;link13;link14;link24;link34;link21;link43];

%% Plot figure
imagesc(networkSpectrum);
colorbar;
xlabel("wavelength channel");
ylabel("link ID");
title('link wavelength blocks');
xlim([1-0.5 8+0.5]);
ylim([1-0.5 7+0.5]);
