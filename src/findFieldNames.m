function Names  = findFieldNames(SimulationSetting, searchPattern)
    Fields  = fieldnames(SimulationSetting);
    nFields = length(Fields);
    
    Names   = cell(nFields, 1);
    nFound  = 0;
    for iField = 1 : nFields
        fName = Fields{iField};
        if regexp(fName, searchPattern, 'once')
            nFound            = nFound + 1;
            Names{nFound}     = fName;
        else
            % skip to next fields
        end
    end;
    Names = Names(1 : nFound);
end
