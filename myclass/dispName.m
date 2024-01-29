function str = dispName(obj)
    keyName  = 'name';
    nObj     = numel(obj);

    strNames = [];
    for iObj = 1 : nObj
        thisObj    = obj(iObj);
        thisString = thisObj.(keyName);
        strNames   = sprintf('%s,%s', strNames, thisString);
    end

    strNames = strNames(2 : end); % skip comma;
    str      = sprintf('%s: [%s]', keyName, strNames);
end
