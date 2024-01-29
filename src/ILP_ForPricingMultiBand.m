function [ConfigNew, obj_PricingProblem,NoLPsConfigNew] = ILP_ForPricingMultiBand( ...
    option_pricing...
    ,ColorlessLightPathSet...
    ,TestBand...
    ,Commodity...
    ,DualVariables ...
    ,Networks...
    )

%=================================
% Parameters to be used;
%=================================
bin_activating_ofp  = logical(TestBand==ColorlessLightPathSet.opticalBandNo);
bin_vec_uv_ofp      = ColorlessLightPathSet.isPathUseEdges';
vec_sefficiency_ofp = ColorlessLightPathSet.capacity;
strPath_ofp         = ColorlessLightPathSet.strPath;
vec_commodityID_ofp = ColorlessLightPathSet.connectionNo(:);


sourceNo_ofconnection       = Commodity.sourceNo;
destinatioNo_ofconnection   = Commodity.destinationNo;
NUM_COMMODITY       = Commodity.nSize;

LinkIndexMatrix     = Networks.edgeIndexMatrix;

MAXIMUM_LIGHTPATHS_CONFIGURE = size(ColorlessLightPathSet.isPathUseEdges,2);

NUM_ROUTES          = length(unique(nonzeros(...
    ColorlessLightPathSet.routeNo...
    )));
NUM_OPTICALBANDS    = length(unique(nonzeros(...
    ColorlessLightPathSet.opticalBandNo...
    )));
NUM_LIGHTPATHS = NUM_COMMODITY * NUM_ROUTES * NUM_OPTICALBANDS;

dual_commodity = DualVariables.('dual_commodity');
dual_onLimitedTransceiver = DualVariables.('dual_onLimitedTransceiver');
dual_onLimitedWavelength  = DualVariables.('dual_onLimitedWavelength');

NoEdge = max(max(LinkIndexMatrix));
%=================================
%          define variable.
%=================================
var_TH          = sdpvar(1,1,'full');
var_bin_ofp     =  binvar(NUM_LIGHTPATHS,1,'full');

%=================================
%          add the constraints.  
%================================ 


Constraints_PricingSubproblem = [];
% Constraints (1) : spectrum non-overlapping constraints;
for lth = 1:NoEdge
    [u,v] = find(LinkIndexMatrix==lth);
    if(~any(bin_vec_uv_ofp(lth,1:NUM_LIGHTPATHS)==1))
        continue;
    end
    Constraints_PricingSubproblem = [Constraints_PricingSubproblem, ...
        (var_bin_ofp(bin_activating_ofp)'*bin_vec_uv_ofp(lth,bin_activating_ofp)'<=1) : sprintf('bandwidthlimit_on(%d,%d)',u,v)];
end


%=================================
%  Set the objective.
%=================================

var_Capacity_oncommodity = sdpvar(NUM_COMMODITY,1,'full');
for idx_commodityi =1: NUM_COMMODITY
    idx_p = bin_activating_ofp & (vec_commodityID_ofp==idx_commodityi);
    var_Capacity_oncommodity(idx_commodityi) = sum(var_bin_ofp(idx_p)'*vec_sefficiency_ofp(idx_p));
end

var_TotalLightpaths   = sum(var_bin_ofp(bin_activating_ofp));

objective_my = (+var_Capacity_oncommodity)'*[dual_commodity']...
    -var_TotalLightpaths*dual_onLimitedTransceiver-dual_onLimitedWavelength(TestBand);

options = sdpsettings('solver', 'gurobi', 'verbose', 0,...
    'cachesolvers',1,...
    'gurobi.FeasibilityTol',1e-9,...
    'gurobi.OptimalityTol',option_pricing.OptimalityTol,...
    'gurobi.MIPGap',option_pricing.MIPgap,...
    'gurobi.TimeLimit', 240);

diagnostics = optimize(Constraints_PricingSubproblem,  -objective_my,options);
 
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
elseif diagnostics.problem == -3
    disp('Solver not found, so of course x is not optimized')
else
    disp('Something else happened')
    disp(diagnostics.info);
end





%=================================
%          Process the data;
%=================================
% ILP - Prcoesse directly the variable;
% % flow_MasterProblem = value(var_flow_p);
bin_ofp                 = value(var_bin_ofp);               % Indicator of LPs;
TotalLightpaths         = value(var_TotalLightpaths);       % Number of LPs 
Capacity_oncommodity    = value(var_Capacity_oncommodity);  % 

obj_PricingProblem      = value(objective_my);

% [model,recoverymodel,diagnostic,internalmodel] = ...
% export(Constraints_PricingSubproblem,objective_my,options);
% display(full([model.A,model.rhs]));
yalmip('clear');
 

% ==============================
% Passing varaibles.
% ==============================

ConfigNew = zeros(MAXIMUM_LIGHTPATHS_CONFIGURE,1);
count = 0;
for idx_commodityi = 1:NUM_COMMODITY
    for p = 1:length(bin_ofp)
        if(bin_ofp(p)==1 && idx_commodityi==vec_commodityID_ofp(p)) % Sorting  the commodity ID.
            count = count + 1;
            ConfigNew(count) = [p];
        end
    end
end

NoLPsConfigNew = TotalLightpaths;

% ==============================
%   Print results.
% ==============================
fprintf('Obj = %d\n',obj_PricingProblem);
fprintf('Capacity_oncommodity = [%s]\n', myPrintArray(Capacity_oncommodity)); 
fprintf('TotalLightpaths = [%d]\n',TotalLightpaths); 

for idx_commodityi = 1:NUM_COMMODITY
    messageCommodity = sprintf('Commodity ID = %d(%d->%d)', ...
        idx_commodityi, ...
        sourceNo_ofconnection(idx_commodityi), ...
        destinatioNo_ofconnection(idx_commodityi)...
        );
    isPrint = false;
    for p = 1:length(bin_ofp)
        if(bin_ofp(p)==1 && idx_commodityi==vec_commodityID_ofp(p))
            cap = vec_sefficiency_ofp(p);
            path = strPath_ofp{p}; 
            messageCommodity = [messageCommodity, ...
                sprintf(', [%s on %d-th band with capacity %g]', ...
                path,...
                TestBand,...
                cap)];
            isPrint = true;
        end
    end
    if(isPrint)
    fprintf('%s\n',messageCommodity);
    else
        ;
    end
end

