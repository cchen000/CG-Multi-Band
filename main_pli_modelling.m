function main_pli_modelling()
% Description:
% 
%  This file presents the SNR modeling accounting for the frequency 
%   dependent loss, nonlinear coefficient, and Raman scattering.
%   - RWBA (MU, ML, MC = -4.4, -4.1, and 0 dB)
%   - RWA  (MU, ML, MC =    0,    0, and 0 dB)
% 
% credited to the following references:
%  [SeKB17] D. Semrau et al. 'Achievable rate degradation of 
%           ultra-wideband coherent fiber communication systems due to 
%           stimulated raman scattering' Opt. Express, vol. 25, no. 12, 
%           pp. 13 024-13 034, Jun. 2017.
%  [ShNS22] N. A. Shevchenko et al. 'Maximizing the information throughput 
%           of ultra-wideband fiber-optic communication systems,' 
%           Opt. Express, vol. 30, no. 11, pp. 19 320-19 331, May 2022
%  [BLBH90] Bredol, M. et al. 'Improved model for OH absorption in optical 
%           fibers.' Journal of Lightwave Technology 8.10 (1990): 1536-1540.
%  [RKHB17] Ian Roberts et al.'Channel Power Optimization of WDM Systems 
%           Following Gaussian Noise Nonlinearity Model in Presence of 
%           Stimulated Raman Scattering' JLT 35.23 (2017): 5237-5249.

% Date: 1st, Jun. 2023.
% Author: cao.chen;

%% SMF fiber setting
LspanSI                     = 100 * 1e3;                            % value : 100 km;
RefLambda                   = 1570 * 1e-9;                          % value : 1550 nm;
LightSpeed                  = 3e8;                                  % value : 30 0000 000 m/s;
RefFrequency                = LightSpeed/RefLambda;                 % m/s;

Btot                        = 15.008 * 1e3 * 1e9;                   % value : 15THz;
Bch                         = 32           * 1e9;                   % value : 32GHz;
NumChannels                 =   Btot/Bch;

NoiseFiguredB               = 5;                                    % value : 5 dB;
NoiseFigure                 = 10^(NoiseFiguredB/10);
PlankConstant               = 6.62607015e-34;       

CenterFrequency_ofch        = [1:NumChannels]'*Bch - Bch/2 - Btot/2; % 1->n, low frequency to high frequency.

% vec_Power0_ofch             = ones(NumChannels,1) * 0.28/4* 1e-3;
vec_Power0_ofch             = ones(NumChannels,1) * 0.448* 1e-3;
Ptot                        = sum(vec_Power0_ofch);

%% Fiber loss
fun_wave                        = @(freq) LightSpeed./(freq);
% Reference [ShNS22]
% Note: with all terms in logarithmic units [dB/km].
% 
RayleighScatteringLossdB        = 1.9 ;                 % [dB/km]. Paper 1.7 dB/km, fix: 1.9 dB/km
RalyeighScatteringWavelength    = 850   * 1e-9;         % [nm]
IRAbsorptionCoefficientdB       = 6.65  * 1e12;         % [dB/km]
IRAbsorptionWavelength          = 52.62 * 1e-6;         % [\mu m];

% ObserveWavelength_wave          = linspace(1400,1700,500)*1e-9;
wave1                           = fun_wave(max(CenterFrequency_ofch)+RefFrequency);
wave2                           = fun_wave(min(CenterFrequency_ofch)+RefFrequency);
ObserveWavelength_wave          = fun_wave(CenterFrequency_ofch +RefFrequency);

Att_Rayleigh                    = RayleighScatteringLossdB * ...
    (RalyeighScatteringWavelength./ObserveWavelength_wave).^4;
Att_IRAbsorption                = IRAbsorptionCoefficientdB*exp(-IRAbsorptionWavelength./ObserveWavelength_wave);

% Peaks (see reference [Improved Model for OH Absorption in Optical Fibers]);
AttPeaks                        = 0.45;
lambda_n                        = [1412e-9      1391e-9     1381e-9     1352e-9     1247e-9];
A_n                             = [0.142        0.606       0.542       0.020       0.059];
sigma_n                     	= [0.00883e6    0.0070e6    0.00440e6   0.00615e6   0.00663e6];

Att3                            = @(lambda)     AttPeaks * A_n(3) * exp( (-1)* ( ( 1./lambda- 1./lambda_n(3))/sigma_n(3)).^2 );
Attn                            = @(n, lambda)  AttPeaks * A_n(n)* 1./( 1+ (( 1./lambda- 1./lambda_n(n))/sigma_n(n)).^2);

Att_Peak                        = Att3(ObserveWavelength_wave)...
    + Attn(1,ObserveWavelength_wave) ...
    + Attn(2,ObserveWavelength_wave) ...
    + Attn(4,ObserveWavelength_wave) ...
    + Attn(5,ObserveWavelength_wave);

Att_wave                        = Att_Rayleigh + Att_IRAbsorption+Att_Peak;

FitAtt0                         = 0.180;
FitAtt1                         = 6.2e-6*(1e9);
FitAtt2                         = 5.7e-6*(1e9)^2;

% Fit for Attenuation;
vec_AttSIFitdB_ofCh               = FitAtt0 ...
    + FitAtt1*(fun_wave(CenterFrequency_ofch+RefFrequency)-fun_wave(RefFrequency)) ...
    + FitAtt2*(fun_wave(CenterFrequency_ofch+RefFrequency)-fun_wave(RefFrequency)).^2;
vec_AttFitSI_ofCh                 = vec_AttSIFitdB_ofCh *log(10)/10/1e3;
vec_AttSI_ofCh                    = vec_AttFitSI_ofCh;

% Fit performance evaluation
figure(1),plot(ObserveWavelength_wave*1e9, [Att_wave,vec_AttSIFitdB_ofCh]);

%% Nonlinear coefficient

% Fit for nonlinear coefficient
FitNonlinearCoefficient0        =   1.2 * 1e-3;
FitNonlinearCoefficient1        =  -2.2 * 1e-6*1e9;

vec_NonlinearCoefficientSI_ofCh =  FitNonlinearCoefficient0 + FitNonlinearCoefficient1*(fun_wave(CenterFrequency_ofch+RefFrequency) - 1570e-9);
%% GVD
GroupVelocityDispersionSignSI   = -21.3 * 1e-3 * (1e-12)^2; % [ps2 km-1]
 
% Raman specified function
bin_FlagISRS                    = 0; % true or false
RamanGainSlopeCoeff             = bin_FlagISRS*0.0290 * 1e-3 * 1e-12; % [1/W/km/THz]
fR                              = 15e3*1e9;               % Raman gain spectrum 15~THz

%% Effective attenuation approach : Raman gain profile derived from [RKHB17, SeKB17]

vec_LeffSI_ofCh                         = (1-exp(-vec_AttSI_ofCh*LspanSI) ) ./vec_AttSI_ofCh;% [RKHB17]

vec_AttEffSRS_ofCh                      = vec_AttSI_ofCh ...
   + 1./LspanSI * (CenterFrequency_ofch * RamanGainSlopeCoeff .* vec_LeffSI_ofCh  .* Ptot ) ...
   + 1./LspanSI * log( fun_sinhc(Ptot/2 * RamanGainSlopeCoeff .* vec_LeffSI_ofCh   * Btot ) );

vec_LeffSRS_ofCh                        = (1-exp(-LspanSI*vec_AttEffSRS_ofCh)) ./vec_AttEffSRS_ofCh;

%% Closed-form approximation evaluation [SeKB17];
% ==============================
% OSNR model in closed-form and AIR calculation;
% ==============================
% OSNR.
NLIEfficiency_ofCh                      = 8/27 * vec_NonlinearCoefficientSI_ofCh.^2 .* vec_AttEffSRS_ofCh .* vec_LeffSRS_ofCh.^2 ...
    .* asinh(0.5 * pi^2 * abs(GroupVelocityDispersionSignSI) * Btot^2 ./vec_AttEffSRS_ofCh ) ...
    ... % Previosuly: asinh(0.5 * pi^2 * abs(GroupVelocityDispersionSignSI) * Btot^2 ./vec_NonlinearCoefficientSI_ofCh ) ...
    ... % Now:        asinh(0.5 * pi^2 * abs(GroupVelocityDispersionSignSI) * Btot^2 ./vec_AttEffSRS_ofCh ) ...
    ./ (pi * abs(GroupVelocityDispersionSignSI)  * Bch^2); 
PNLIClosedForm_ofCh                     = NLIEfficiency_ofCh .*vec_Power0_ofch.^3;
PASEClosedForm_ofCh                     = NoiseFigure * PlankConstant  * RefFrequency * Bch .* exp(vec_AttEffSRS_ofCh * LspanSI);

OpticalSNR0ClosedForm_ofCh              = vec_Power0_ofch./( PASEClosedForm_ofCh + PNLIClosedForm_ofCh);

figure(2),clf;
hold on;
plot(CenterFrequency_ofch * 1e-12,10*log10(OpticalSNR0ClosedForm_ofCh));
xlim([min(CenterFrequency_ofch)* 1e-12 max(CenterFrequency_ofch)* 1e-12]);
xlabel('Channel Frequency [THz], f_w');
ylabel('SNR [dB]');

%% Margin computation.
BandID                          = zeros(NumChannels,1);
B1                              = floor(NumChannels/3);
B2                              = floor(NumChannels/3);
B3                              = NumChannels - B1 - B2;
BandID(1:B1)                    = 1;
BandID((B1+1):(B1+B2))          = 2;
BandID((B2+B1+1):(B1+B2+B3))    = 3;

WorstCaseSNR                    =  min(OpticalSNR0ClosedForm_ofCh);
WorstCaseSNRB1                  =  min(OpticalSNR0ClosedForm_ofCh(BandID==1));
WorstCaseSNRB2                  =  min(OpticalSNR0ClosedForm_ofCh(BandID==2));
WorstCaseSNRB3                  =  min(OpticalSNR0ClosedForm_ofCh(BandID==3));

fprintf('Power=%g dBm/ch, WorstCase SNR = %.3g dB, Mb=(%.2g,%.2g,%.2g)dB\n', ...
    10*log10(vec_Power0_ofch(1)) + 30,...
    10*log10(WorstCaseSNR), ...
    10*log10(WorstCaseSNR/WorstCaseSNRB1), ...
    10*log10(WorstCaseSNR/WorstCaseSNRB2), ...
    10*log10(WorstCaseSNR/WorstCaseSNRB3));

end
function y = fun_sinhc(x)
y = (eps^3+sinh(x))./(eps^3+x);

end