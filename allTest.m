%%
%  *test gurobi*
addpath(genpath('src\'));
addpath(genpath('myclass\'));
addpath(genpath('testcase\'));
yalmiptest();

%%
%  *test user-function*
a1 = runtests('testcase\TransceiverTest.m');
a2 = runtests('testcase\ConfigurationTest.m');
a3 = runtests('testcase\getTest.m');
a4 = runtests('testcase\printTest.m');

nFailed = sum([a1.Failed, a2.Failed, a3.Failed, a4.Failed]);

%%
assert(nFailed == 0 ...
    , 'MATLAB:test:Error' ...
    , 'test failed');
