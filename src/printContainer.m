function [varargout] = printContainer(mapObj)
    %export keywords and keyvalues of map container.
    % Input:
    %   - container
    % Output(varargout):
    %   - nameRow
    %   - valueRow
    %
    %   Example:
    %         keySet =   {'Jan', 'Feb', 'Mar', 'Apr'};
    %         valueSet = [327.2, 368.2, 197.6, 178.4];
    %         dataRow = containers.Map(keySet,valueSet);
    %         printContainer(dataRow);
    %
    % ==============================
    %   Obtain values in container
    % ==============================
    itemNames  = mapObj.keys;
    itemValues = mapObj.values;
    nItem      = mapObj.Count;
    % ==============================
    %   Export data
    % ==============================
    nameRow    = [];
    for iNameItem  = 1 : nItem
        cache   = sprintf('%s;', itemNames{iNameItem});
        nameRow = sprintf('%s%s', nameRow, cache);
    end
    
    valueRow  = [];
    for iValueItem = 1 : nItem
        if ischar(itemValues{iValueItem})
            cache    = sprintf('%s;',  itemValues{iValueItem});
            valueRow = sprintf('%s%s', valueRow, cache);
        elseif isnumeric(itemValues{iValueItem})
            cache    = sprintf('%g;',  itemValues{iValueItem});
            valueRow = sprintf('%s%s', valueRow, cache);
        else
            error('not defined');
        end
    end
    
    % ==============================
    % Deal with function output
    % ==============================
    nOutputs  = nargout;
    varargout = cell(1, nOutputs);
    
    for iOutput = 1 : nOutputs
        switch iOutput
            case 1
                varargout{1} = nameRow;
            case 2
                varargout{2} = valueRow;
            otherwise
                error('undefined');
        end
    end

end