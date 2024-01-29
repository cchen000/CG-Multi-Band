classdef ConfigurationTest < matlab.unittest.TestCase
    methods(Test)
        
        
        function sortConfigTest1(testCase)
            
            ConfigNo    = [1 2 3 4 5 6];
            repeatConfigurations = [1 0 2 0 0 3];
            sort_metric = [1 0 2 0 0 3]; % repeatConfigurations
            
            % actSolution
            [configNoThisBand, repConfigNo, repTimes] = sortConfig(ConfigNo ...
                , repeatConfigurations...
                , sort_metric);

            % expSolution
            expRepConfigNo = [6 3 1];
            expRepTimes    = [3 2 1];
            expConfigNoThisBand = [6 6 6 3 3 1];

            % verify
            testCase.verifyEqual(expRepConfigNo, repConfigNo);
            testCase.verifyEqual(expRepTimes, repTimes);
            testCase.verifyEqual(expConfigNoThisBand,configNoThisBand);
        end
        
        function sortConfigTest2(testCase)
            ConfigNo    = [1 2 3 4 5 6];
            repeatConfigurations = [1 0 2 0 0 3];
            sort_metric = [20 1 10 100 100 1000];
            
            % actSolution
            [configNoThisBand, repConfigNo, repTimes] = sortConfig(ConfigNo ...
                , repeatConfigurations...
                , sort_metric);

            % expSolution
            expRepConfigNo = [6 1 3];
            expRepTimes    = [3 1 2];
            expConfigNoThisBand = [6 6 6 1 3 3];

            % verify
            testCase.verifyEqual(expRepConfigNo, repConfigNo);
            testCase.verifyEqual(expRepTimes, repTimes);
            testCase.verifyEqual(expConfigNoThisBand,configNoThisBand);
        end
        
        function sortConfigTest3(testCase)
            
            ConfigNo    = [1 2 3 4 5 6];
            repeatConfigurations = [1 0 2 0 0 3];
            sort_metric = [20 1 1000 100 100 10];
            
            % actSolution
            [configNoThisBand, repConfigNo, repTimes] = sortConfig(ConfigNo ...
                , repeatConfigurations...
                , sort_metric);

            % expSolution
            expRepConfigNo = [3 1 6];
            expRepTimes    = [2 1 3];
            expConfigNoThisBand = [3 3 1 6 6 6];

            % verify
            testCase.verifyEqual(expRepConfigNo, repConfigNo);
            testCase.verifyEqual(expRepTimes, repTimes);
            testCase.verifyEqual(expConfigNoThisBand,configNoThisBand);
        end
        
        
        function sortConfigTest4(testCase)
            
            ConfigNo    = [1 6 5 3 2 4];
            repeatConfigurations = [1 3 0 2 0 0];
            sort_metric = [20 1 1000 100 100 10];
      
            % actSolution
            [configNoThisBand, repConfigNo, repTimes] = sortConfig(ConfigNo ...
                , repeatConfigurations...
                , sort_metric);
            
            % expSolution
            expRepConfigNo = [3 1 6];
            expRepTimes    = [2 1 3];
            expConfigNoThisBand = [3 3 1 6 6 6];

            % verify
            testCase.verifyEqual(expRepConfigNo, repConfigNo);
            testCase.verifyEqual(expRepTimes, repTimes);
            testCase.verifyEqual(expConfigNoThisBand,configNoThisBand);
        end
        
    end
end

