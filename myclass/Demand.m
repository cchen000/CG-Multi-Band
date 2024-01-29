classdef Demand < matlab.mixin.Copyable
% Demand: This class denotes a set that may contain multiple entities.
%     'nSize'           : [1] size of sequential demands, |SEQ|, i\in SEQ.
%     'no'              : [matrix |SEQ|*1] ID of the i-th demand
%     'nodePairNo'      : [matrix |SEQ|*1] nodePairNo of i-th demand
%     'sourceNo'        : [matrix |SEQ|*1] source node of i-th demand
%     'destinationNo'   : [matrix |SEQ|*1] destination node of i-th demand
%     'flowSize'        : [matrix |SEQ|*1] traffic required by i-th demand
    properties
        nSize           = [];
        no              = [];
        nodePairNo      = [];
        sourceNo        = []; 
        destinationNo   = []; 
        flowSize        = [];
    end

    methods 
        
        function obj = Demand(nElements)
            narginchk(1,1);
            if(~isnumeric(nElements))
                error('Input must be numerical');
            end
            
            obj.nSize           = 0;
            obj.no              = zeros(nElements,1,'uint64');
            obj.nodePairNo      = zeros(nElements,1,'uint16');
            obj.sourceNo        = zeros(nElements,1,'uint8');
            obj.destinationNo   = zeros(nElements,1,'uint8');
            obj.flowSize        = zeros(nElements,1,'double');
        end
        
        function [cTargetNo] = findNodePairNo(obj, vertexList)
            % find the number of vertexList matching demand <obj>
            s = vertexList(1);
            d = vertexList(2);
            cTargetNo = [];
            for iNodePair = 1 : obj.nSize
                if (obj.sourceNo(iNodePair)==s)&&...
                        (obj.destinationNo(iNodePair)==d)
                    cTargetNo = obj.nodePairNo(iNodePair);
                    break;
                end
            end
        end
        
        function obj_plus = appendn(obj1, obj2, n)
            narginchk(3,3);
            if(~isa(obj1,'Demand') || ...
                    ~isa(obj2,'Demand'))
                error('Input must be demand types');
            end
            if(obj2.nSize >n)
                error('Size does not match!');
            end
            nBasicSize = obj1.nSize;
            obj_plus = copy(obj1);
            
            keyName         ='nSize';
            obj_plus.(keyName) = obj1.(keyName) + n;
            
            PropertyLists   = properties(obj1);
            
            for i = 1:length(PropertyLists)
                props = PropertyLists{i};
                if(strcmp(props, keyName))
                    continue;
                end
                
                DimSize =  size(obj2.(props));
                if(DimSize(2)<=1&&length(DimSize)<=2)
                    obj_plus.(props)((nBasicSize+1):(nBasicSize+n)) = ...
                        obj2.(props)(1:n);
                else
                   error('Check syntax for non-column copy'); 
                end
            end
        end
        
        function obj2 = indexing(obj,idx)
            narginchk(2,2);
            KeyNameSize        = 'nSize';
            nIDX               = length(idx);
            obj2               = Demand(length(idx));
            obj2.(KeyNameSize) = length(idx);
            
            PropertiesList      = properties(obj);
            
            for i = 1:length(PropertiesList)
                props = PropertiesList{i};
                if(strcmp(props,KeyNameSize))
                    continue;
                end
                obj2.(props)(1:nIDX) = obj.(props)(idx(1:nIDX));
            end
        end
    end
end
