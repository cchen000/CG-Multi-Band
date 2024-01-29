
function [obj_RestrictedMasterProblem, Repeated_Configurations] = solveRMP_INT( ...
    OptionRMPINT,...
    CandidateWaveConfigSet, ...
    Commodity, vec_limitedwavelengths_onbands,LimitedTransceivers...
    )
% same as the solveRMP.m except that varaibles are integer;
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
var_Repeated_Configurations    =  intvar(NoConfigs,1,'full');

%=================================
%          add the constraints.  
%================================ 

Constraints_MaxCapWDM = [];

% Constraints (1) : Network throughput definition.
for idx_commodityi = 1:NUM_COMMODITY
    Constraints_MaxCapWDM = [Constraints_MaxCapWDM, ...
        ( ... 
        commodity_amount(idx_commodityi)/sum(commodity_amount) * var_TH  ...
        + (...
        - Capacity_commodityi_onconfigc(idx_commodityi,bin_activating_ofconfig) ...
        * [var_Repeated_Configurations(bin_activating_ofconfig)])...
        <= 0) : sprintf('commodity_%d',idx_commodityi)];
end

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

options = sdpsettings('solver', 'gurobi', 'verbose', 2,...
    'cachesolvers',1,...
    'gurobi.FeasibilityTol',1e-9,...
    'gurobi.OptimalityTol',OptionRMPINT.OptimalityTol,...
    'gurobi.MIPGap',OptionRMPINT.MIPgap,...
    'gurobi.TimeLimit', OptionRMPINT.TimeLimit);

diagnostics = optimize(Constraints_MaxCapWDM, - objective_my,options);
 
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
% 
%          Process the data;
% 
%=================================
% ILP - Prcoesse directly the variable;
% flow_MasterProblem = value(var_flow_p);
Repeated_Configurations          = value([var_Repeated_Configurations]);
Repeated_Configurations(isnan(Repeated_Configurations)) = 0;
obj_RestrictedMasterProblem      = value(objective_my);

Configurations_onconfig = groupSum(...
    Repeated_Configurations(:), ...
    vec_bandtypes_onconfigc);

vec_LPs_onoband = groupSum(...
    NoLPs_ofconfig(:).*Repeated_Configurations(:), ...
    vec_bandtypes_onconfigc);

% dual_commodity         = zeros(1,NUM_COMMODITY);
% for idx_commodity = 1:NUM_COMMODITY
%     dual_commodity(idx_commodity) = dual(Constraints_MaxCapWDM(sprintf('commodity_%d',idx_commodity)));
% end
% dual_onLimitedTransceiver = dual(Constraints_MaxCapWDM('Limited LPs'));
% dual_onLimitedWavelenghth  = dual(Constraints_MaxCapWDM('Limited Wavelengths'));

% [model,recoverymodel,diagnostic,internalmodel] = ...
% export(Constraints_MaxCapWDM,objective_my,options);
% display(full([model.A,model.rhs]));
yalmip('clear');

% fprintf('**Restricted Master Problem==============================\n');
fprintf('Obj = %d\n',obj_RestrictedMasterProblem);
fprintf('No. Times = [%s]\n',printArray(Repeated_Configurations(:)'));
fprintf('opticalBandNo = [%s]\n', printArray(vec_bandtypes_onconfigc));
% fprintf('Dual of commodities : [%s]\n',num2str(dual_commodity));
% fprintf('Dual of limited transceivers : [%s]\n',num2str(dual_onLimitedTransceiver));
% fprintf('Dual of limited wavelength   : [%s]\n',num2str(dual_onLimitedWavelenghth));
fprintf('No. Transceivers = sum([%s]) = [%s]/[%s]\n',printArray(vec_LPs_onoband), num2str(NoLPs_ofconfig * [Repeated_Configurations]), num2str(LimitedTransceivers)); 
% fprintf('No. Wavelengths  = [%s]/[%s]\n',num2str(sum(1* [Repeated_Configurations])), num2str(NUM_CHANNELS_perBands')); 
fprintf('No. Wavelengths  = [%s]/[%s]\n',printArray(Configurations_onconfig), printArray(vec_limitedwavelengths_onbands)); 
