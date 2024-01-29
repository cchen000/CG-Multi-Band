function str = printCell(obj, delimeter)
    narginchk(1,2);
    
    if nargin == 1
        delimeter = ',';
    else
        % deal with two inputs
    end
    
    % 'cell'
    assert(isa(obj, 'cell') ...
        , 'MATLAB:printCell:InvalidTypes' ...
        , 'Wrong input'...
        );
    
    str  = [];
    nObj = length(obj);
    for iObj = 1 : nObj
        thisObj = obj{iObj};
        if isa(thisObj, 'char')
            cache   = sprintf('%s', thisObj);
            str     = sprintf('%s%s%s', str, cache, delimeter);
        elseif isa(thisObj, 'double')
            cache   = sprintf('%g', thisObj);
            str     = sprintf('%s%s%s', str, cache, delimeter);
        else
            error('undefined');
        end
    end
    
    lenStr                                         = length(str);
    str((lenStr - length(delimeter) + 1) : lenStr) = [];
end