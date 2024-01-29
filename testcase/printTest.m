classdef printTest < matlab.unittest.TestCase
    methods(Test)
        function printContainerTest(testCase)
            % actSolution;
            keySet =   {'Jan', 'Feb', 'Mar', 'Apr'};
            valueSet = [327.2, 368.2, 197.6, 178.4];
            dataRow = containers.Map(keySet,valueSet);
            [actNameRow, actValueRow] = printContainer(dataRow);
            
            % expSolution;
            expNameRow  = 'Apr;Feb;Jan;Mar;';
            expValueRow = '178.4;368.2;327.2;197.6;';
            
            % verify
            testCase.verifyEqual(actNameRow,expNameRow)
            testCase.verifyEqual(actValueRow,expValueRow)
        end
        function printArrayTest(testCase)
            % actSolution;
            valueSet = [327.2, 368.2, 197.6, 178.4];
            actSolution = printArray(valueSet);
            
            % expSolution
            expSolution = '327.2,368.2,197.6,178.4';
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
        
        function printCellTest(testCase)
            % actSolution;
            valueSet    = {327.2, 368.2, 197.6, 178.4};
            actSolution = printCell(valueSet);
            
            % expSolution
            expSolution = '327.2,368.2,197.6,178.4';
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
    end
end

