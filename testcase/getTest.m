classdef getTest < matlab.unittest.TestCase
% test functions with prefix 'get*';
% get* function can access the contents of a class obj.
    methods(Test)
        
        function getNumeric1Test(testCase)
            % structure test
            s = struct(...
                'Jan',      327.2 ...
                , 'Feb',    368.2 ...
                , 'Mar',    197.6 ...
                , 'Apr',    178.4 ...
                );
            % actSolution;
            actSolution = getNumeric(s, 'Jan');
            
            % expSolution;
            expSolution = 327.2;
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
        
        
        function getNumeric2Test(testCase)
            % class test
            
            % actSolution;
            actSolution = getNumeric(testclass(), 'Jan');
            
            % expSolution;
            expSolution = 327.2;
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
        
        
        function getChar1Test(testCase)
            % class test
            
            % actSolution;
            actSolution = getChar(testclass(), 'Zodiac');
            
            % expSolution;
            expSolution = 'snake';
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
        
        
        
        function getChar2Test(testCase)
            % structure test
            s = struct(...
                'Jan',      327.2 ...
                , 'Feb',    368.2 ...
                , 'Mar',    197.6 ...
                , 'Apr',    178.4 ...
                , 'Zodiac', 'snake'...
                );
            % actSolution;
            actSolution = getChar(s, 'Zodiac');
            
            % expSolution;
            expSolution = 'snake';
            
            % verify
            testCase.verifyEqual(actSolution,expSolution)
        end
        
        
        function getClassByFields1Test(testCase)
            names = 'Smith, Daniel, Jane';
            keyName = 'name';
            totalObj(1) = struct(...
                'name', 'Smith'...
                , 'employNo', 1 ...
                , 'age', 20 ...
                );
            totalObj(2) = struct(...
                'name', 'Daniel'...
                , 'employNo', 2 ...
                , 'age', 23 ...
                );
            totalObj(3) = struct(...
                'name', 'Jack'...
                , 'employNo', 3 ...
                , 'age', 28 ...
                );
            totalObj(4) = struct(...
                'name', 'Jane'...
                , 'employNo', 4 ...
                , 'age', 25 ...
                );
            
            [candidateObj, nCandidate] = getClassByFields(names, keyName, totalObj);
            
            
            testCase.verifyEqual(nCandidate,3)
            testCase.verifyEqual(candidateObj(1),totalObj(1));
            testCase.verifyEqual(candidateObj(2),totalObj(2));
            testCase.verifyEqual(candidateObj(3),totalObj(4));
        end
    end
end

