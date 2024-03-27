classdef WaveConfig < matlab.mixin.Copyable
%Wavelength Configuration set that contains multiple configurations;
    properties
        nSize                 = 0; 
        no                    = []; % format : 1 x |nConfigs|
        isActive              = []; % format : 1 x |nConfigs|
        nLightPath            = []; % format : 1 x |nConfigs|
        nRepetition           = []; % format : 1 x |nConfigs|
        nFrequencySlot        = []; % format : 1 x |nConfigs|
        opticalBandNo         = []; % format : 1 x |nConfigs|
        lightPathNoSet        = []; % format : |MAX_LIGHTPATHS| x |nConfigs|
        capacity_perCommodity = []; % format : |nCommodity| x |nConfigs|
    end
    
    methods
        
        function str = disp(obj)
            FieldsName = fieldnames(obj);
            
            str = [];
            for iName = 1:numel(FieldsName)
                thisName = FieldsName{iName};
                if(isa(obj.(thisName),'char'))
                    str = sprintf('%s%s[%s]: %s\n',str, thisName,printArray(size(val), '*'), obj.(keyName));
                elseif(isnumeric(obj.(thisName))&&numel(obj.(thisName))==1)
                    val = obj.(thisName);
                    str = sprintf('%s%s[%s]: [%s]\n',str, thisName,printArray(size(val), '*'), sprintf("%g",val));
                elseif(isnumeric(obj.(thisName))&&numel(obj.(thisName))>=2&&size(obj.(thisName),1)==1)
                    val = obj.(thisName);
                    str = sprintf('%s%s[%s]: [%s]\n',str, thisName,printArray(size(val), '*'), printArray(val));
                elseif(isnumeric(obj.(thisName)) &&numel(obj.(thisName))>=2&&size(obj.(thisName),1)>=2)
                    val = obj.(thisName);
                    [x1,y1] = find(val);
                    resultStr = [];
                    nElement  = min(40, length(x1));
                    for iElement = 1:nElement
                        cacheStr = printArray(val(x1(iElement), y1(iElement)));
                        resultStr = sprintf('%s{%s}', resultStr, cacheStr);
                    end
                    str = sprintf('%s%s[%s](%d elements): [%s]\n',str, thisName,printArray(size(val), '*'), nElement, resultStr);
                    % do not deal with full matrix;
                elseif(isa(obj.(thisName),'logical')&&size(obj.(thisName),1)==1)
                    val = obj.(thisName);
                    str = sprintf('%s%s[%s]: [%s]\n',str, thisName, printArray(size(val), '*'), printArray(val));
                elseif(isa(obj.(thisName),'logical') && size(obj.(thisName),1)>=2)
                    val = obj.(thisName);
                    [x1,y1] = find(val);
                    resultStr = [];
                    nElement  = min(40, length(x1));
                    for iElement = 1:nElement
                        cacheStr = printArray([x1(iElement), y1(iElement)]);
                        resultStr = sprintf('%s{%s}', resultStr, cacheStr);
                    end
                    str = sprintf('%s%s[%s](%d elements): [%s]\n',str, thisName,printArray(size(val), '*'), nElement, resultStr);
                    % do not deal with full matrix;
                elseif(isa(obj.(thisName),'cell'))
                    val = obj.(thisName);
                    resultStr = [];
                    nElement  = min(20, length(x1));
                    for iElement = 1:nElement
                        cacheStr = sprintf('%s',val{iElement});
                        resultStr = sprintf('%s{%s}', resultStr, cacheStr);
                    end
                    str = sprintf('%s%s[%s](%d elements): %s\n',str, thisName,printArray(size(val), '*'),nElement, resultStr);
                else
                    error('not defined yet');
                end
            end
        end
        
        function obj2 = ConfigureParameters(obj, ...
                configNoSet, ...
                LP_Config_Block,...
                ColorlessLightPathSet)
            % Description:
            %
            %       Configure parameters based on lightpath no. set;
            %
            % ==============================
            % Import parameters to be used later;
            % ==============================
            obj2            = copy(obj);
            if(isempty(LP_Config_Block))
                return; 
            end
            nConfigs        = length(configNoSet);
            assert(nConfigs<=size(LP_Config_Block,2), 'Wrong size, please expand a larger space');
            NUM_COMMODITY  = length(unique(nonzeros(...
                ColorlessLightPathSet.nodePairNo...
                )));
            vec_hops_ofp        = ColorlessLightPathSet.hops;
            vec_capacity_ofp    = ColorlessLightPathSet.capacity;
            vec_commodityNo_ofp = ColorlessLightPathSet.nodePairNo;
            vec_bandNo_ofp = ColorlessLightPathSet.opticalBandNo;
            % ==============================
            % Assign parameters
            % ==============================
            
            % set Active;
            for iConfig = 1:nConfigs
                cNo = configNoSet(iConfig);
                obj2.isActive(cNo) = true;
                obj2.no(cNo) = cNo;
            end
            
            % Assign size;
            obj2.nSize = nnz(obj2.no);
            
            % Assign pNo and related parameters;
            % :: give p_{s,d,k,b} ID to c
            % :: compute a_{c}, ob_{c} based on p_{s,d,k,b};
            for iConfig = 1:nConfigs
                cNo = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                assert(...
                    length(pNoSet)==length(unique(pNoSet)),...
                    'Repeated LPs');
                
                nSet = nnz(pNoSet);
                obj2.lightPathNoSet(1:nSet,cNo) = pNoSet;
            end
            for iConfig = 1:nConfigs
                cNo = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                
                obj2.nLightPath(cNo) = nnz(pNoSet);
            end
            for iConfig = 1:nConfigs
                cNo = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                
                obNo = unique(vec_bandNo_ofp(pNoSet));
                assert(length(obNo)==1,'Wrong:: not in the same band');
                obj2.opticalBandNo(cNo) = obNo;
            end
            
            
            
            % Assign frequency slots;
            for iConfig   = 1:nConfigs
                cNo    = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                
                obj2.nFrequencySlot(cNo) = sum(vec_hops_ofp(pNoSet));
            end
            
            
            
            % Assign capacity;
            % :: compute T_{s,d,k,c}(= sum_{k} delta_{s,d,k,c}* C_{s,d,k},\forall s,d,c
            for iConfig = 1:nConfigs
                cNo    = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                [cap,targetConnectionNo] = groupSum(...
                    vec_capacity_ofp(pNoSet)...
                    , vec_commodityNo_ofp(pNoSet) ...
                    );
                obj2.capacity_perCommodity(targetConnectionNo,cNo) = cap;
            end
            
            
            % ==============================
            % Check wavelength clash
            % ==============================
            for iConfig = 1:nConfigs
                cNo    = configNoSet(iConfig);
                pNoSet = nonzeros(LP_Config_Block(:,iConfig));
                
                assert(all(sum(ColorlessLightPathSet.isPathUseEdges(pNoSet,:),1)<=1), ...
                    'Error: Wavelength clash!');
            end
            
        end
        

        function obj = WaveConfig(nConfigs,paras)
            narginchk(2,2);
            
            assert(isa(paras,'struct'),'Not proper input');
            
            MAX_POSSIBILITY_PERCONFIGURE = paras.MAXIMUM_LIGHTPATHS_CONFIGURE;
            NUM_COMMODITY = paras.NUM_COMMODITY;
            assert(...
                isnumeric(nConfigs) ...
                &&isnumeric(MAX_POSSIBILITY_PERCONFIGURE)...
                &&isnumeric(NUM_COMMODITY), 'Value must be numeric');
            
            obj.nSize = nConfigs;
            obj.isActive = false(1,nConfigs);
            obj.nLightPath = zeros(1,nConfigs,'uint16');
            obj.nRepetition = zeros(1,nConfigs,'uint64');
            obj.nFrequencySlot = zeros(1,nConfigs);
            obj.lightPathNoSet = zeros(MAX_POSSIBILITY_PERCONFIGURE, nConfigs,'uint64');
            obj.capacity_perCommodity = zeros(NUM_COMMODITY,nConfigs, 'double');
            obj.opticalBandNo = zeros(1, nConfigs, 'uint8');
        end
        
        function [obj] = updateConfigurationTimes(obj ...
                , rep_Config)
            obj.nRepetition = rep_Config;
        end
        
        function [vec_ColoredConfigID_ofCh] = packConfiguration(...
                  CandidateWaveConfigSet ...
                , PackOption ...
                )
            % Description:
            %       This function assigns the exact wavelength for each configuration.
            %
            % Default: Repeat each configuration based on configraution times.
            %
            % Date: 20 Nov. 2023.
            % Author: cao.chen
            % ==============================
            % Obtain parameters to be used later;
            % ==============================
            isGroup = PackOption.('isGroup');
            direction = PackOption.('direction');
            basis  = PackOption.('basis');            
            Repeated_Configurations = CandidateWaveConfigSet.nRepetition;
            nTotalTimes = sum(Repeated_Configurations);
            nConfigs = CandidateWaveConfigSet.nSize;

            random_perconfig = randperm(nConfigs);
            % Random swapping before processing,
            %    in order to ensure different sorting results;
            % Or, simply set as 1:nConfigs;
            
            NUM_BANDS = length(unique(nonzeros(...
                CandidateWaveConfigSet.opticalBandNo(random_perconfig)...
                )));
            nFS_perconfig = CandidateWaveConfigSet.nFrequencySlot(random_perconfig);
            nLPs_perconfig = CandidateWaveConfigSet.nLightPath(random_perconfig);
            vec_bandtypes_perconfig = CandidateWaveConfigSet.opticalBandNo(random_perconfig);

            Repeated_Configurations = Repeated_Configurations(random_perconfig);
            
            
            vec_ColoredConfigID_ofCh = [];
            % ==============================
            % Pack and repeat configurations on each optical band;
            % ==============================
            for ob = 1:NUM_BANDS
                
                % find the activating configurations;
                bin_config = ...
                    (Repeated_Configurations(:)>=1)...
                    &(vec_bandtypes_perconfig(:)==ob);
                
                % choose sorting metrics;
                switch(basis)
                    case 'numberOfRepetitions'
                        sort_metric      = Repeated_Configurations;
                    case 'numberOfFrequencySlots'
                        sort_metric      = nFS_perconfig;
                    case 'numberOfLightPaths'
                        sort_metric      = nLPs_perconfig;
                    case 'random'
                        sort_metric      = randperm(nConfigs);
                    otherwise
                        error('non-specified sorting options');
                end
                
                % sort and expand configuration by configTimes;
                [configNoThisBand,~,~] = ...
                    sortConfig(random_perconfig, Repeated_Configurations(:).*bin_config(:),sort_metric, ...
                    direction);
                
                vec_ColoredConfigID_ofCh = [vec_ColoredConfigID_ofCh, configNoThisBand];
            end
            
            if(isGroup==true)
                ; % do nothing as previous is already group sorted;
            elseif(isGroup==false && strcmp(basis,'random'))
                NewChannelOrder   = randperm(nTotalTimes);
                vec_ColoredConfigID_ofCh = vec_ColoredConfigID_ofCh(NewChannelOrder);
            else
                error('not defined');
            end
        end
        
        function printConfiguration(    ...
                obj                     ...
                , fileName              ...
                , itemConfiguration     ...
                , itemLightPath         ...
                , ColorlessLightPathSet)
            % Description:
            %       prints wavelength configuration to a file.
            %
            % Date: 30 Nov. 2023.
            % Author: cao.chen
            % - Print information from the given items, 26th Mar., 2024.
            % ==============================
            % Obtain parameters
            % ==============================
            NUM_COMMODITY                = size(obj.capacity_perCommodity,1);
            MAXIMUM_LIGHTPATHS_CONFIGURE = size(ColorlessLightPathSet.isPathUseEdges,2);
            rep                          = obj.nRepetition; % nRepetition_ofConfig
            
            % Start print ...
            fid = fopen(fileName,'a+');
            fprintf(fid,'==============================\n');
            
            count = 0;
            for iConfig = 1 : numel(rep)
                if rep(iConfig) == 0
                    continue; 
                end 
                count   = count + 1;
                % ==============================
                % print configurations ...
                % ==============================
                fprintf(fid, '%d-th configuration (', count);
                for jItemConfig = 1 : numel(itemConfiguration)
                    thisItem = itemConfiguration{jItemConfig};
                    val      = obj.(thisItem);
                    
                    fprintf(fid, '%s=[%s],' ...
                        , thisItem ...
                        , printArray(val(:,iConfig))...
                        );
                end 
                fprintf(fid, ')\n');
                 
                % ==============================
                % and print lightpaths ...
                % ==============================
                pNoSet           = nonzeros(...
                    obj.lightPathNoSet(1:MAXIMUM_LIGHTPATHS_CONFIGURE,iConfig)...
                    ); 
                 
                thisLightPathSet = ...
                    LightPath(...
                    MAXIMUM_LIGHTPATHS_CONFIGURE ...
                    , size(ColorlessLightPathSet.isPathUseEdges,2) ...
                    ); 
                for jPath = 1 : length(pNoSet)
                    pNo              = pNoSet(jPath);
                    thisLightPathSet = ...
                        addLightPath(...
                          thisLightPathSet ...
                        , ColorlessLightPathSet ...
                        , pNo);
                end 
                 
                printLightPathSet(...
                    thisLightPathSet ...
                    , itemLightPath ...
                    , fileName, 'a+'...
                    )
            end
            
            fclose(fid);
            
        end
    end
end