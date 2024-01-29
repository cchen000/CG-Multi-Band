function str = getChar(obj, keyName)
    % get properties of <keyName> in array class <obj>
    narginchk(1,2);
    if nargin == 1
        keyName = 'name';
    end
    assert(any(strcmp(fieldnames(obj(1)), keyName))...
        , 'MATLAB:getValue:InvalidProperty' ...
        , 'No such properties');
    assert(ischar(obj(1).(keyName)) ...
        , 'MATLAB:getValue:TypeError'...
        , 'This function only supports chars');
    str = printCell(...
        arrayfun(@(thisObj) thisObj.(keyName) ...
        , obj ...
        , 'UniformOutput',false) ...
        );
end