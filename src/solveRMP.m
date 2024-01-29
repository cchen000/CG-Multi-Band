function [obj_RestrictedMasterProblem, Repeated_Configurations, dual_commodity, dual_onLimitedTransceiver, dual_onLimitedWavelength,ReducedCost_perConfig] ...
    = solveRMP( ...
    CandidateWaveConfigSet, ...
    Commodity, vec_limitedwavelengths_onbands,LimitedTransceivers,...
    OptionsForRMP)
% Description:
%       Restricted Master Problem;  
% <CandidateWaveConfigSet> : candidate wavelength configuration set;
% <Commodity> : network demand informations;
% vec_limitedwavelengths_onbands: max. wavelengths of optical band;
% LimitedTransceivers:  max. transceivers of network;
% ==============================
% Obtain parameters to be used later;
% ==============================
commodity_amount               = Commodity.flowSize;
NUM_COMMODITY                  = Commodity.nSize;
NoConfigs                      = CandidateWaveConfigSet.nSize;
Capacity_commodityi_onconfigc  = CandidateWaveConfigSet.capacity_perCommodity(...
    1:NUM_COMMODITY,1:NoConfigs);
NoLPs_ofconfig                 = double(CandidateWaveConfigSet.nLightPath(1:NoConfigs));
vec_bandtypes_onconfigc        = CandidateWaveConfigSet.opticalBandNo(1:NoConfigs);
bin_activating_ofconfig        = CandidateWaveConfigSet.isActive(1:NoConfigs);
NUMBER_OF_BANDS                = length(vec_limitedwavelengths_onbands);
% ==============================
%  Define varaiables.
% ==============================

var_TH = sdpvar(1,1,'full');
var_Repeated_Configurations    =  sdpvar(NoConfigs,1,'full');

%=================================
%          add the constraints.  
%================================ 

Constraints_MaxCapWDM = [];

% Constraints (1) : Network throughput definition.
% for idx_commodityi = 1:NUM_COMMODITY
%     Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
%         ( ... 
%         commodity_amount(idx_commodityi)/sum(commodity_amount) * var_TH  ...
%         + (-Capacity_commodityi_onconfigc(idx_commodityi,bin_activating_ofconfig)*[var_Repeated_Configurations(bin_activating_ofconfig)]) ...
%         <=0 ) : sprintf('commodity_%d',idx_commodityi)];
% end
Constraints_NetworkThroughput = [( ...
    commodity_amount(:) /sum(commodity_amount) * var_TH ...
    +  (-Capacity_commodityi_onconfigc(1:NUM_COMMODITY,bin_activating_ofconfig) *[var_Repeated_Configurations(bin_activating_ofconfig)])<=0 )];

% Constraints (2) : limited LPs.
Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
    (NoLPs_ofconfig(bin_activating_ofconfig) * [var_Repeated_Configurations(bin_activating_ofconfig)] ...
    <= LimitedTransceivers) : 'Limited LPs'];

% % Constraints (3) : limited wavelengths.
% Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
%     (sum([var_Repeated_Configurations])<=LimitedWavelength):'Limited Wavelengths'];

for b = 1:NUMBER_OF_BANDS
    idx_bandtypes = (vec_bandtypes_onconfigc==b);
    if(any(idx_bandtypes))
        Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
            (sum([var_Repeated_Configurations(idx_bandtypes & bin_activating_ofconfig)])<=vec_limitedwavelengths_onbands(b)):sprintf('Limited Wavelengths on %d-th band',b)];
    end
end


% Constraints (4) : positive value;
Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
    ([var_Repeated_Configurations(bin_activating_ofconfig)]>=0)];

%=================================
%          Set the objectives.  
%================================ 

objective_my = full(var_TH);


diagnostics = optimize([Constraints_MaxCapWDM;Constraints_NetworkThroughput], - objective_my,OptionsForRMP);
 
if diagnostics.problem == 0
%     disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
elseif diagnostics.problem == -3
    disp('Solver not found, so of course x is not optimized')
else
    disp('Something else happened')
    disp(diagnostics.info);
end




%=================================
% 
%          Process the data;
% 
%=================================
% ILP - Prcoesse directly the variable;
% flow_MasterProblem = value(var_flow_p);
Repeated_Configurations          = value([var_Repeated_Configurations]);
obj_RestrictedMasterProblem      = value(objective_my);

% ===Backup Code for variables <dual_commodity>===
% dual_commodity    = zeros(1,NUM_COMMODITY);
% for idx_commodity = 1:NUM_COMMODITY
%     dual_commodity(idx_commodity) = dual(Constraints_MaxCapWDM(sprintf('commodity_%d',idx_commodity)));
% end
dual_commodity = dual(Constraints_NetworkThroughput)';

dual_onLimitedTransceiver = dual(Constraints_MaxCapWDM('Limited LPs'));
dual_onLimitedWavelength  = zeros(NUMBER_OF_BANDS,1);
for b = 1:NUMBER_OF_BANDS 
    try
        dual_onLimitedWavelength(b) = dual(Constraints_MaxCapWDM(sprintf('Limited Wavelengths on %d-th band',b)));
    catch
        dual_onLimitedWavelength(b) = 0;
    end
end

% Reduced Cost for existing configurations.
ReducedCost_perConfig = -1e5 + zeros(1,NoConfigs);
for i = 1:NoConfigs
    b = vec_bandtypes_onconfigc(i);
    ReducedCost_perConfig(i) = dual_commodity * Capacity_commodityi_onconfigc(1:NUM_COMMODITY,i) - dual_onLimitedWavelength(b)  - NoLPs_ofconfig(i) * dual_onLimitedTransceiver;
end

% [model,recoverymodel,diagnostic,internalmodel] = ...
% export(Constraints_MaxCapWDM,objective_my,options);
% display(full([model.A,model.rhs]));
yalmip('clear');

% fprintf('**Restricted Master Problem==============================\n');
% fprintf('Obj = %.15e\n',obj_RestrictedMasterProblem);
% fprintf('Repeated Times = [%s]\n',num2str(Repeated_Configurations')); 
% fprintf('Conf. ID = [%s]\n', num2str(vec_bandtypes_onconfigc));
% fprintf('(Dual Constr.) ReducedCost = [%s]\n',num2str(ReducedCost_perConfig)); 
% fprintf('Dual of commodities : [%s]\n',num2str(dual_commodity));
% fprintf('Dual of limited transceivers : [%s]\n',num2str(dual_onLimitedTransceiver));
% fprintf('Dual of limited wavelength   : [%s]\n',num2str(dual_onLimitedWavelength'));
