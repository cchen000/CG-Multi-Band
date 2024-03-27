classdef LightPath < matlab.mixin.Copyable
%denotes a lightpath set containing multiple lightpaths;
%       
% LightPath 
%   'nSize'         : [1] number of lightpaths, |P|, p_{s,d,k,w}\in P.
%   'no'            : [integer |P|*|E|] ID of a lightpath p_{s,d,k,w}
%   'nodePairNo'    : [integer |P|*1] n_(s,d), node pair number of p_{s,d,k,w}
% 
%   'sourceNo'      : [integer |P|*1] s for path p_{s,d,k,w};
%   'destinationNo' : [integer |P|*1] d for path p_{s,d,k,w};
%   'routeNo'       : [integer |P|*1] k for path p_{s,d,k,w}.
%   'isPathUseEdges': [integer |P|*|E|] links required by path p_{s,d,k,w}.
% 
%   'strPath'       : [cell string |P|*1] links string of path;
%   'cost'          : [matrix |P|*1] length of a lightpath;
%   'hops'          : [integer |P|*1] number hops of a lightpath;
% 
%   'transmissionModeNo': [integer |P|*1] transmission mode no. of the ...
%                     transceviers in a lightpath;
%   'transceiverNo' : [integer |P|*1] transmission mode no. of the ...
%                     transceviers in a lightpath;
%   'nFrequencySlot': [integer |P|*1] no. slots of a lightpath
%   'wavelengthNo'  : [integer |P|*1] wavelength no. of a lightpath
%   'opticalBandNo' : [integer |P|*1] optical band no.of a lightpath
% 
%   'capacity' : [|P|*1] transmission capacity of a lightpath
% 
    properties
        nSize       = 0;
        nodePairNo  = [];
        no          = [];
        
        sourceNo        = [];
        destinationNo   = [];
        routeNo         = [];
        isPathUseEdges  = [];

        strPath         = {};
        cost            = [];
        hops            = [];
        
        %  Parameters (trans, modulation formats, oband, etc. );
        transmissionModeNo = [];
        transceiverNo      = []; 
        nFrequencySlots    = [];
        wavelengthNo       = [];
        opticalBandNo      = [];
        
        capacity = [];
    end
    
    methods
        
        function disp(obj)
            if(isempty(obj))
                str = 'containing empty LightPath\n';
                fprintf('%s\n',str);
            elseif(numel(obj)==1)
                str = dispElement(obj);
                fprintf('%s\n',str);
            elseif(numel(obj)>=2)
                str = dispArray(obj);
                fprintf('%s\n',str);
            else
                str = 'not defined';
                fprintf('%s\n',str);
            end
        end
        
        function str = dispElement(obj)
            % display detailed information for a obj.
            fNames = fieldnames(obj);
            str = [];
            for iName = 1:numel(fNames)
                thisName = fNames{iName};
                if(isa(obj.(thisName),'char'))
                    str = sprintf('%s%s: %s\n',str, thisName, obj.(thisName));
                elseif(isnumeric(obj.(thisName))&&numel(obj.(thisName))==1)
                    val = obj.(thisName);
                    str = sprintf('%s%s: [%s]\n',str, thisName,sprintf("%g",val));
                elseif(isnumeric(obj.(thisName))&&numel(obj.(thisName))>=2)
                    val = obj.(thisName);
                    str = sprintf('%s%s: [%s]\n',str, thisName,printArray(val));
                elseif(isa(obj.(thisName),'logical') && length(size(obj.(thisName)))==2)
                    val = obj.(thisName);
                    [x1,y1] = find(val);
                    resultStr = [];
                    nElement = min(40,length(x1));
                    for iElement = 1:nElement
                        cacheStr = printArray([x1(iElement), y1(iElement)]);
                        resultStr = sprintf('%s{%s}', resultStr, cacheStr);
                    end
                    str = sprintf('%s%s(%d elements): [%s]\n',str, thisName,nElement, resultStr);
                    % do not deal with mat;
                elseif(isa(obj.(thisName),'cell'))
                    val = obj.(thisName);
                    resultStr = [];
                    nElement = min(20, length(val));
                    for iElement = 1:nElement
                        cacheStr = sprintf('%s',val{iElement});
                        resultStr = sprintf('%s{%s}', resultStr, cacheStr);
                    end
                    str = sprintf('%s%s(%d elements): %s\n',str, thisName, nElement, resultStr);
                elseif(isempty(obj.(thisName)))
                    str = sprintf('%s%s: %s\n',str, thisName, 'empty');
                else
                    error('not defined yet');
                end
            end
        end
        
        function flag = isempty(obj)
            % true if all elements are empty,
            % false otherwise,
            n       = numel(obj);
            flag    = true;
            for i   = 1:n
                if(obj(i).nSize~=0)
                    flag = false;
                end
            end
        end
        function str = dispArray(obj)
            % output size information for array objects.
            keyName = 'nSize';
            str  = sprintf('%s of these sets:',keyName);
            for i = 1:numel(obj)
                isPrint =false;
                if(obj(i).(keyName)>=1)
                    isPrint = true;
                end
                if(isPrint)
                    str = sprintf('%s{%d:%d}',str,i,obj(i).(keyName));
                end
            end
        end
        
        
        function [varargout] = convertLightPath2EdgeSlot(...
                LightPathSet, isChannelUseSlots)
            %
            % Input:
            %	isChannelUseSlots          indicator between channel and slot
            %	LightPathSet               lightpath set
            % Output(varargout):
            %	edgeSlotBlock              edge slot status
            %	edgeSlotBlock_withPathNo   edge slot status, by PathNo
            %
            % ==============================
            % Obtain parameters to be used
            % ==============================
            pNoSet          = nonzeros(LightPathSet.no);
            isPathUseEdges  = LightPathSet.isPathUseEdges(pNoSet, 1 : end);
            channelNoOfp    = LightPathSet.wavelengthNo(pNoSet);
            
            % ==============================
            % define basic parameters
            % ==============================
            nSlots          = size(isChannelUseSlots, 2);
            isPathUseSlots  = isChannelUseSlots(channelNoOfp, 1 : nSlots);
            
            % ==============================
            % Output selective results
            % ==============================
            nOutputs        = nargout;
            varargout       = cell(1, nOutputs);
            
            for iOutput = 1 : nOutputs
                switch iOutput
                    case 1
                        edgeSlotBlock               = double(isPathUseEdges)'...
                            * double(isPathUseSlots);
                        varargout{1}                = edgeSlotBlock;
                    case 2
                        edgeSlotBlock_withPathNo    = double(isPathUseEdges)'...
                            * diag(pNoSet) ...
                            * double(isPathUseSlots);
                        varargout{2}                = edgeSlotBlock_withPathNo;
                    otherwise
                        error('undefined');
                end
            end
        end

        function printLightPathSet(self ...
                , itemName ...
                , fileName ...
                , FormatSpec)
            %       prints LightPathSet to a file.
            NUM_SET       = numel(self);
            
            fid = fopen(fileName,FormatSpec);
            
            % short outline
            if NUM_SET >= 2
            fprintf(fid, '%d LightPathSets, total LPs =%d\n'...
                , NUM_SET ...
                , sum(getNumeric(self, 'nSize')) ...
                );
            end
            
            for iSet = 1:NUM_SET
                thisSet = self(iSet);
                
                % detailed information;
                for jPath = 1 : thisSet.nSize
                    fprintf(fid, ':::');
                    cache           = '';
                    for kItem = 1 : length(itemName)
                        assert(...
                            ismember(itemName{kItem}, fieldnames(thisSet))...
                            , 'MATLAB:LightPath:UnMatchedName'...
                            , 'error name'...
                            );
                        
                        thisName    = itemName{kItem};
                        thisItem    = thisSet.(thisName);
                        if iscell(thisItem)
                            thisVal = thisItem{jPath};
                        elseif isnumeric(thisItem)
                            thisVal = printArray(thisItem(jPath));
                        else
                            error('undefined');
                        end
                        cache = sprintf('%s %s=%s,',cache, thisName, thisVal);
                    end
                    fprintf(fid , '%s\n', cache);
                end
            end
            fclose(fid);
        end
        
        function printLightPathSetWithFlow(self ...
                , Commodity ...
                , itemName ...
                , fileName ...
                , FormatSpec)
            % Description:
            %
            %       prints LightPathSet to a file.
            %
            % ==============================
            % Date: 30 Nov. 2023.
            % Author: cao.chen
            NUM_SET       = numel(self);
            
            fid = fopen(fileName,FormatSpec);
            
            % short outline
            fprintf(fid, '%d LightPathSets, total LPs =%d\n'...
                , NUM_SET ...
                , sum(getNumeric(self, 'nSize')) ...
                );
            
            for iSet = 1:NUM_SET
                thisSet = self(iSet);
                
                cNo     = unique(nonzeros(thisSet.nodePairNo)); 
                assert(...
                    ismember(cNo,Commodity.no)&&length(cNo)==1,...
                    'Wrong commodity index');
                
                % summary information;
                nLPs	= thisSet.nSize;
                fprintf(fid,...
                    '%d-th commodity (%d-->%d, flow=%.3g, no. LPs=%d, cap=%g)\n', ...
                    cNo, ...
                    Commodity.sourceNo(cNo), ...
                    Commodity.destinationNo(cNo), ...
                    Commodity.flowSize(cNo), ...
                    thisSet.nSize, ...
                    sum(thisSet.capacity(1:nLPs))...
                    );
                                
                % detailed information;
                for jPath = 1 : thisSet.nSize
                    fprintf(fid, ':::');
                    cache           = '';
                    for kItem = 1 : length(itemName)
                        assert(...
                            ismember(itemName{kItem}, fieldnames(thisSet))...
                            , 'MATLAB:LightPath:UnMatchedName'...
                            , 'error name'...
                            );
                        
                        thisName    = itemName{kItem};
                        thisItem    = thisSet.(thisName);
                        if iscell(thisItem)
                            thisVal = thisItem{jPath};
                        elseif isnumeric(thisItem)
                            thisVal = printArray(thisItem(jPath));
                        else
                            error('undefined');
                        end
                        cache = sprintf('%s %s=%s,',cache, thisName, thisVal);
                    end
                    fprintf(fid , '%s\n', cache);
                end
            end
            fclose(fid);
        end
            
        function obj = LightPath(val,nEdges)
            
            if(nargin==0)
                return; % allow empty object initilization;
            end
            
            narginchk(2,2);
            if(isnumeric(val)&&isnumeric(nEdges))
                obj.nodePairNo    = zeros(val, 1, 'double');
                obj.no            = zeros(val, 1, 'double');
                
                obj.sourceNo      = zeros(val, 1, 'uint16');
                obj.destinationNo = zeros(val, 1, 'uint16');
                obj.routeNo       = zeros(val, 1, 'uint8');
                obj.isPathUseEdges= zeros(val, nEdges, 'logical');
                
                obj.hops          = zeros(val, 1, 'uint16');
                obj.cost          = zeros(val, 1, 'uint16');
                obj.strPath       = cell(val, 1);
                
                obj.transmissionModeNo = zeros(val, 1, 'uint8');
                obj.transceiverNo      = zeros(val, 1, 'uint8'); 
                obj.nFrequencySlots    = zeros(val, 1, 'double');
                obj.wavelengthNo       = zeros(val, 1, 'uint64');
                obj.opticalBandNo      = zeros(val, 1, 'uint8');
                
                obj.capacity           = zeros(val, 1, 'double');
            else
                error('Value must be numeric')
            end
        end
        
        function obj = addLightPath(obj1, obj2, pNoSet)
            % Description:
            %      Add lightpaths pNoSet from obj2 to obj1.
            % ==============================
            % Date: 30th, Nov. 2023.
            % Author: cao chen;
            
            % Check arguments;
            narginchk(3,3);
            assert(isa(obj1,'LightPath')&&isa(obj2,'LightPath'), ...
                'not proper LightPath');
            
            % Check size;
            sizeName     = 'nSize';
            nSizeBasic   = obj1.nSize;
            nSizeApp     = numel(pNoSet);
            assert(obj2.(sizeName)>=nSizeApp,'Not supported appendix');
            
            % Make a copy from obj1 to obj, 
            %   then move elements from obj2 to obj..
            obj          = copy(obj1);
            fNames       = properties(obj1);
            for i = 1:length(fNames)
                thisName = fNames{i};
                
                if strcmp(thisName, sizeName)
                    obj.(sizeName) = obj1.(sizeName) + nSizeApp;
                    continue;
                else 
                    % copy other properties
                end
                
                dimSize =  size(obj2.(thisName));
                if(dimSize(2)==1||dimSize(1)==1) % default: copy all;
                    obj.(thisName)((nSizeBasic+1):(nSizeBasic+nSizeApp)) = ...
                        obj2.(thisName)(pNoSet); % not efficient for large n.
                elseif(strcmp(thisName, 'isPathUseEdges'))
                    [~, NoEdges] = size(obj1.isPathUseEdges);
                    obj.isPathUseEdges((nSizeBasic+1):(nSizeBasic + nSizeApp) ,1:NoEdges)...
                        = obj2.isPathUseEdges(pNoSet,1:NoEdges);
                else
                    error('Check syntax for non-column copy');
                end
            end
        end
        
        function obj2 = createColorlessLightPath(...
                obj ...
                , NodePair ...
                , NetState ...
                )
            %create a colorless lightpath set, used within a configuraion
            % Note:
            %   - path is considered as p_{s,d,k,b};  
            %   - default assumption: transceiverNo is 1;
            % 
            % ==============================
            % Date: 18 Nov. 2023
            % Author: cao.chen
            obj2 = copy(obj);
            % ==============================
            % pass parameter
            % ==============================
            bin_sdk_atl   = NetState.isPathOnsdkUseLink;
            candidatePath = NetState.candidatePathOfsdk;
            costOfCandidatePath = NetState.lengthOfPathOfsdk;
            nOpticalBands = max(NetState.opticalBandNoOfSlot);
            
            NUM_NODEPAIR = NodePair.nSize;
            [NUM_NODE,~,NUM_ROUTE, NUM_EDGE] = size(bin_sdk_atl);
            assert(NUM_NODE*(NUM_NODE-1)>=NUM_NODEPAIR,'Wrong commodity');
            % ==============================
            % Computation starts ...
            %       default transceiverNo is 1;
            % ==============================
            pNo = 0;
            for iNodePair = 1:NUM_NODEPAIR
                s = NodePair.sourceNo(iNodePair);
                d = NodePair.destinationNo(iNodePair);
                for k = 1:NUM_ROUTE
                    if(any(bin_sdk_atl(s,d,k,1:NUM_EDGE)==1))
                        
                        for ob = 1:nOpticalBands
                            pNo             = pNo + 1;
                            obj2.no(pNo)    = pNo;
                            obj2.nodePairNo(pNo) = iNodePair;
                            
                            obj2.sourceNo(pNo)      = s;
                            obj2.destinationNo(pNo) = d;
                            obj2.routeNo(pNo)       = k;
                            obj2.isPathUseEdges(pNo, 1 : NUM_EDGE) ...
                                = bin_sdk_atl(s, d, k, 1 : NUM_EDGE);

                            obj2.strPath{pNo}   = ...
                                printArray(candidatePath{s, d, k}, '-');
                            obj2.cost(pNo) ...
                                                = costOfCandidatePath(s, d, k);
                            obj2.hops(pNo)      = ...
                                sum(bin_sdk_atl(s, d, k, 1 : NUM_EDGE));
                            
                            obj2.opticalBandNo(pNo)= ob;
                            obj2.transceiverNo(pNo) = 1; % 1 default
                        end
                    end
                end
            end
            
            obj2.nSize = pNo; % number starts from 1, rather than 0;            
            obj2.no((pNo+1):end) = [];
            obj2.sourceNo((pNo+1):end) = [];
            obj2.destinationNo( (pNo+1):end) = [];
            obj2.routeNo((pNo+1):end) = [];
            obj2.nodePairNo((pNo+1):end) = []; 
            obj2.isPathUseEdges((pNo+1):end,:) = [];
            
            obj2.opticalBandNo((pNo+1):end) = [];
            obj2.strPath((pNo+1):end)= [];
            obj2.hops((pNo+1):end) = [];
            
            obj2.transceiverNo((pNo+1):end) = [];
        end
        
    end
end