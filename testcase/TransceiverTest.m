classdef TransceiverTest < matlab.unittest.TestCase
    methods(Test)
        function initialTest(testCase)
            valueArray = [10, 10, 1, 20, -2100, 10];
            ValueStructure = struct('name', 'Trans10G'...
                ,'bandwidth', 10 ...
                ,'capacity', 20 ...
                ,'numberOfFrequencySlots', 1 ...
                ,'penalty', -2100 ...
                ,'baud_rate', 10 ...
                );
            
            % actSolution
            act_Trans  = Transceiver(ValueStructure);
            
            % expSolution
            exp_Trans  = Transceiver();
            exp_Trans.name = 'Trans10G';
            exp_Trans.bandwidth = 10;
            exp_Trans.nFrequencySlots = 1;
            exp_Trans.GBaud = 10;
            exp_Trans.penalty = -2100;
            exp_Trans.exampleCapacity = 20;
            
            % verify
            testCase.verifyEqual(act_Trans,exp_Trans)
        end
    end
end

