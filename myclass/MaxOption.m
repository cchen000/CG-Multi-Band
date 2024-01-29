classdef MaxOption
    % Show properties value by closing function disp();
    properties
        algorithmName = '',
        algorithmOption = struct([]),
        rngSeed = 0;
        costFcn = {};
        costLimit = 1e5; % Ideally, much would be better for infinity.
        capacityScalingFactor = 1;
    end
    methods
        
        function obj2 = setSolverOption(obj,varargin)
            nFields = (nargin-1)/2;
            assert(rem(nFields,1)==0,'Input should be a pair');
            propertyTable = reshape(varargin,[2,nFields])';
            for iField = 1:nFields
                fieldname = strsplit(propertyTable{iField,1},".");
                value     = propertyTable{iField,2};
                if(numel(fieldname)==1)
                    assert(isprop(obj,fieldname{1}),'Not defined property name');
                    obj.(fieldname{1}) = value;
                elseif(numel(fieldname)==2)
                    assert(isprop(obj,fieldname{1}),'Not defined property name');
                    assert(isfield(obj.(fieldname{1}),fieldname{2}),'Not defined field name');
                    obj.(fieldname{1}).(fieldname{2}) = value;
                else
                    error('not defined');
                end
            end
            obj2 = obj;
        end
        
        function obj = MaxOption(solverName)
            assert(isa(solverName,'char'),'Input must be char');
            defaultILPOption = struct(...
                'algorithmName', 'ILPmodel', ...
                'ILP_FeasibilityTol', 1e-6, ... SimulationSetting.ILP_FeasibilityTol, ...
                'ILP_OptimalityTol', 1e-6, ... SimulationSetting.ILP_OptimalityTol, ...
                'ILP_MIPGap', 2e-2, ... SimulationSetting.ILP_mingap, ...
                'ILP_TimeLimit', 3600); %SimulationSetting.ILP_TimeLimit
            defaultColumnGenerationOption = struct(...
                'algorithmName', 'ColumnGenerationHeu', ...
                'strategyLPConfigInitial', '1LP per configuration', ...
                ...%'1LP per configuration', 'kSP-FF per configuration'
                'CG_NoLoops', 400, ... SimulationSetting.iteration_CG, ...
                'rmp_mingap', 1e-2, ... SimulationSetting.rmp_mingap, ...
                'rmp_TimeLimit', 10, ... SimulationSetting.rmp_TimeLimit, ...
                'RMPINTOptimalityTol', 1e-2, ...  SimulationSetting.rmpINT_optimality_tol, ...
                'RMPINTMIPgap', 1e-2, ... SimulationSetting.rmpINT_mingap, ...
                'RMPINTTimeLimit', 10, ... SimulationSetting.rmpINT_TimeLimit,...
                'ConfigGroupSorting', 'number_of_Reps',...
                ...%  'number_of_FSs', 'number_of_LPs', 'number_of_Reps', 'randomly',
                'ConfigGroupDirection', 'descend', ... %  'ascend' or 'descend',
                'ConfigSingleDirection', 'default',...%  'randomly', 'default'
                'pricing_optimality_tol', 1e-4, ... SimulationSetting.pricing_optimality_tol, ...
                'pricing_mingap', 5e-2); % SimulationSetting.pricing_mingap
            defaultSequentialOption  = struct( ...
                'algorithmName', 'SequentialLoadingHeu', ...
                'LoadingStrategy','kSP-FF', ...
                'TerminateStrategy', 'First-blocking');
            
            % Set an algorithm option based on <solverName>
            switch(solverName)
                case 'SequentialLoadingHeu'
                    obj.algorithmOption=defaultSequentialOption;
                case 'ColumnGenerationHeu'
                    obj.algorithmOption=defaultColumnGenerationOption;
                case 'ILPmodel'
                    obj.algorithmOption=defaultILPOption;
                otherwise
                    error('Not defined solver');
            end
            obj.algorithmName = solverName;
        end
    end
end