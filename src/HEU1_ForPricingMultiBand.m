function [ConfigNew, obj_PricingProblem,TotalLightpaths] = HEU1_ForPricingMultiBand( ...
    option_pricing, ColorlessLightPathSet,TestBand,Commodity,DualVariables)
% Maximum Weight Independent Set Problem, Weight = | Reduced Cost | .
% :: Heuristic algorithm 1, select the node based on the weight's descending order.

bin_activating_ofp  = logical(TestBand==ColorlessLightPathSet.opticalBandNo);
%=================================
%       Read parameters;
%=================================
NUM_ROUTES          = uint64(max(ColorlessLightPathSet.routeNo));
NUM_OPTICALBANDS    = uint64(max(ColorlessLightPathSet.opticalBandNo));
MAXIMUM_LIGHTPATHS_CONFIGURE = size(ColorlessLightPathSet.isPathUseEdges,2);
NoEdge              = size(ColorlessLightPathSet.isPathUseEdges,2);

vec_commodityID_ofp = ColorlessLightPathSet.nodePairNo;
capacity_ofp        = ColorlessLightPathSet.capacity;
strPath_ofp         = ColorlessLightPathSet.strPath;
bin_vec_uv_ofp      = ColorlessLightPathSet.isPathUseEdges';

NUM_COMMODITY               = uint64(Commodity.nSize);
sourceNo_ofconnection       = Commodity.sourceNo;
destinatioNo_ofconnection   = Commodity.destinationNo;

NUM_LIGHTPATHS = NUM_COMMODITY * NUM_ROUTES * NUM_OPTICALBANDS;
% Note: store as uint8 is fine, 
%      but needs conversion when assessing variables.

dual_commodity = DualVariables.('dual_commodity');
dual_onLimitedTransceiver = DualVariables.('dual_onLimitedTransceiver');
dual_onLimitedWavelength  = DualVariables.('dual_onLimitedWavelength');

%=================================
%          Parameter Definition
%=================================
ReducedCost_ofp     = zeros(NUM_LIGHTPATHS,1);
bin_vec_ofp         = zeros(NUM_LIGHTPATHS,1);    % Indicator of LPs;
bin_vec_uvofc       = zeros(NoEdge,1);
Capacity_oncommodity= zeros(NUM_COMMODITY,1);

% Calculate ReducedCost/lightpath
for pid = 1:NUM_LIGHTPATHS
    if(bin_activating_ofp(pid)==0) % Only considers the activated lightpath.
        continue;
    end
    commodityID = vec_commodityID_ofp(pid);
    if(commodityID>0)
        ReducedCost_ofp(pid) = dual_commodity(commodityID) * ( + capacity_ofp(pid)) ...
            - dual_onLimitedTransceiver;
    end
end

% Sort lightpath p by its Reduced Cost.
[~,ascendIndex_ofp] = sort(ReducedCost_ofp,'descend');

% Calculate ReducedCost of Pricing Problem.
ReducedCostAccum = 0;
for idx_p   = 1:NUM_LIGHTPATHS
    pid                     = ascendIndex_ofp(idx_p);
    if(bin_activating_ofp(pid)==0) ... % Only considers the activated lightpath.
        continue;
    end
    if(ReducedCost_ofp(pid)>=1e-4)
    bin_uv                  =  bin_vec_uv_ofp(1:NoEdge, pid);
    % Select pid if ...
    if(max(bin_vec_uvofc + bin_uv)<=1 ...  % No Wavelength-Clash
            && 1)         % Negative Reduced Cost.
        bin_vec_uvofc       = bin_vec_uvofc + bin_uv;
        bin_vec_ofp(pid)	= true;
        ReducedCostAccum    = ReducedCostAccum ...
            + ReducedCost_ofp(pid);
    end
    else
        break;
    end
end
obj_PricingProblem = ReducedCostAccum + (-dual_onLimitedWavelength(TestBand));

% ==============================
% Post-process 
% ==============================
% :: Calculate Total LPs.
% :: Calculate Respective Capacity on a Commodity.
TotalLightpaths    = sum(bin_vec_ofp);
for pid = 1:NUM_LIGHTPATHS
    % :: Respective capacity on a connection.
    if(bin_vec_ofp(pid)==true)
        commodityID = vec_commodityID_ofp(pid);
        Capacity_oncommodity(commodityID) = Capacity_oncommodity(commodityID) ...
            + capacity_ofp(pid);
    end
end

% ==============================
% Passing varaibles.
% ==============================

ConfigNew = zeros(MAXIMUM_LIGHTPATHS_CONFIGURE,1);
count = 0;

LPs_IDset   = find(bin_vec_ofp);
for idx_commodityi = 1:NUM_COMMODITY
    % Sorting lightpathID based on commodity ID.
    for idx_p = 1:length(LPs_IDset)
        pid     = LPs_IDset(idx_p);
        if(idx_commodityi==vec_commodityID_ofp(pid)) % Sorting  the commodity ID.
            count = count + 1; 
            ConfigNew(count) = pid;
        end
    end
end

% ==============================
%   Print results.
% ==============================
% fprintf('Obj = %d\n',obj_PricingProblem);
% fprintf('Capacity_oncommodity = [%s]\n', myPrintArray(Capacity_oncommodity)); 
% fprintf('TotalLightpaths = [%d]\n',TotalLightpaths); 
% 
% for idx_commodityi = 1:NUM_COMMODITY
%     messageCommodity = sprintf('Commodity ID = %d(%d->%d), ', ...
%         idx_commodityi, ...
%         sourceNo_ofconnection(idx_commodityi), ...
%         destinatioNo_ofconnection(idx_commodityi)...
%         );
%     isPrint = false;
%     for p = 1:length(bin_vec_ofp)
%         if(bin_vec_ofp(p)==1 && idx_commodityi==vec_commodityID_ofp(p))
%             cap = capacity_ofp(p);
%             path = strPath_ofp{p}; 
%             messageCommodity = [messageCommodity, ...
%                 sprintf('[%s on %d-th band with capacity %g],', ...
%                 path,...
%                 TestBand,...
%                 cap)];
%             isPrint = true;
%         end
%     end
%     if(isPrint)
%     fprintf('%s\n',messageCommodity);
%     else
%         ;
%     end
% end

