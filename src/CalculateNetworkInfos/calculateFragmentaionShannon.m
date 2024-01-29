function H = calculateFragmentaionShannon(A)
% 
% Title: Minimum- and Maximum-Entropy Routing and Spectrum Assignment for Flexgrid Elastic Optical Networking
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
%         T =[1 1 1 1 1 0 1 0 0 1 0 0 0 1 1 1];
%         A = T;
% 
%  Created by Cao CHEN (chen.cao{at}sjtu.edu.cn)
%  Date: 16-09-2021

D = length(A);
i = 0;
H = 0;


Astr = num2str(A);
Astrblock = split(Astr,{'1'});
for i = 1:length(Astrblock)
    NonZeros = find(double(Astrblock{i,:})~=' ');
    if(~isempty(NonZeros))
        Di = length(NonZeros);
        H = H + Di/D * log(Di/D);
    end
end

H = - H;


