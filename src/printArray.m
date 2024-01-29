function str = printArray(val, delimeter)
    
    if nargin == 1
        delimeter = ',';
    else
        % deal with two inputs
    end

    expectedFormat = sprintf('%%g%s', delimeter);
    str            = sprintf(expectedFormat, val);
    
    lenStr                                         = length(str);
    str((lenStr - length(delimeter) + 1) : lenStr) = [];
end