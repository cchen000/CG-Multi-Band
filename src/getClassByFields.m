function [candidateObj, nCandidate] = getClassByFields(names, keyName, totalObj)
%get candidates of total, according to keyword <keyName>.
    candidateName            = splitstring(names);
    nCandidate               = length(candidateName);
    candidateObj(nCandidate) = totalObj(1);    
    assert(...
        length(totalObj) >= nCandidate...
        , 'MATLAB:getClassByFields:OutOfRange'...
        , 'Candidate should be less than total'...
        );
    for iCandidate = 1 : nCandidate
        nFound = 0;
        for jObj = 1 : length(totalObj)
            thisObj = totalObj(jObj);
            if strcmp(thisObj.(keyName), candidateName{iCandidate})
                candidateObj(iCandidate) = thisObj;
                nFound                   = nFound + 1;
            else
                % skip to next candidate, as this obj does not coincide
            end
        end;
        
        assert(nFound == 1 ...
            , 'MATLAB:getClassByFields:WrongFound'...
            , 'none existing or duplicate parameters');
    end
end

function itemName = splitstring(names)
    names    = replace(names, " ", ""); % remove whitespace
    itemName = strsplit(names, {','});
end
