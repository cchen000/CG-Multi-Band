function H = calculateFragmentaionABP(A,G)
% 
% Title: Spectrum fragmentation issue in flexible optical networks: analysis and good practices
% Description: 
%       This function calculates the Entropy of a fiber link based on the
%       reference [1].
% 
% Example: 
%         A = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%         A1 =[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
%         A2 =[1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0];
%         A3 =[1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0];
%         A4 =[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
%         A5 =[1 1 1 1 0 0 1 0 0 1 0 0 0 1 1 1];
% 
%         S1 =[1 1 1 0 0 0 1 0 0 1 0 0 0 1 1 1];
%         S2 =[1 1 0 0 0 0 1 0 0 1 0 0 0 0 1 1];
%         S3 =[1 1 0 0 0 0 1 0 0 1 0 0 0 0 1 1];
%         S4 =[1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0];
% 
%         T = [1 1  0 0 1 1 0 0 1 1];
%         G = [2];
%         A = T;
% 
%  Created by Cao CHEN (chen.cao{at}sjtu.edu.cn)
%  Date: 21-09-2021

Astr = num2str(A);
Astrblock = split(Astr,{'1'});

% -Calculate Numerator
Numerator = 0; 
for i = 1:length(Astrblock)
    NonZeros = find(double(Astrblock{i,:})~=' ');
    if(~isempty(NonZeros))
        Di = length(NonZeros);
        for k = 1:length(G)
            Numerator = Numerator + floor(Di/G(k));
        end
    end
end

% -Calculate Denominator 
Denominator = 0; 
Di = 0;
for i = 1:length(Astrblock)
    NonZeros = find(double(Astrblock{i,:})~=' ');
    if(~isempty(NonZeros))
        Di = Di + length(NonZeros);
    end
end

for k = 1:length(G)
    Denominator = Denominator + floor( Di/G(k) );
end


if(Numerator==0 && Denominator==0) % In case of full occupation.
    H = 0;
else
H = 1 - Numerator / Denominator;
end


