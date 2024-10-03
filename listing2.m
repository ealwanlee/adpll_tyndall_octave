%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PLL MODEL WHICH INCLUDES %%%%
%%% A TDC WITH QUANTIZATION %%%%
%%% AND DITHER %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters obtained from %%%%%
%%%%% PLL Assistant Designer %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fz/fo = 0.22
% order 2
% shape Butter
% Ref freq 10e6
% out freq = 1e9
% fz = 2.2e5;
% wz = 2*pi*fz;
% fp = 1.734e5;
% wp = 2*pi*fp;
% N = 100;
% Kccpsim = 7.272e10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nsamples = 2^22;
Nsamples = 2e5;
fz = 2.2e4;
wz = 2*pi*fz;
fp = 1.734e5;
wp = 2*pi*fp;
fref = 10e6;
N = 100;
Kccpsim = 7.272e10;
KccpsimAdapted = Kccpsim/wz;
kv = 10e6;
kp = 1;
klf = N*KccpsimAdapted/(kv*kp);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inizialization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
freq_freerunning = fref*N;
output_phase_dco = zeros(Nsamples ,1);
phase_input = zeros( Nsamples ,1);
phase_err = zeros( Nsamples ,1);
freq_dco = zeros( Nsamples ,1);
out = zeros(Nsamples ,1);
outExtra = zeros( Nsamples ,1);
outExtraGain = zeros(Nsamples ,1);
input_dco = zeros( Nsamples ,1);
output_phase_dco (1) = 1;
freq_dco(1) = freq_freerunning;
eq = zeros(1,Nsamples);
tin = zeros(1,Nsamples);
outTDC = zeros(1,Nsamples);
outTDCwithGain = zeros(1, Nsamples);
eqDeltaC = zeros(1, Nsamples);
Finput = fref * ones(Nsamples ,1); % Input Frequency
time_array =(1/fref) * (0:1: Nsamples -1); % Time Array
fs = fref; % Sampling Frequency of the Digital Filter
ditherDCO = zeros(1,length (time_array));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN : ZERO plus INTEGRATOR %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z z1 outz inz p1;
s = 2 * fs * (z-1)/(z+1);
Hz = (s + z1) / (s);
Hz = simplify(Hz);
[N,D] = numden(Hz);
Cnum = coeffs(N, z);
Cden = coeffs(D, z);
Cnum1 = Cnum(1);
Cnum2 = Cnum(2);
Cden1 = Cden(1);
Cden2 = Cden(2);
fsval = 2*pi*fs; % sampling frequency in the discrete time domain
z1val = wz;
Cnum1val = double( subs(Cnum1 , {z1 , fs},{z1val , fsval}));
Cnum2val = double( subs(Cnum2 , {z1 , fs},{z1val , fsval}));
Cden1val = double( subs(Cden1 , {z1 , fs},{z1val , fsval}));
Cden2val = double( subs(Cden2 , {z1 , fs},{z1val , fsval}));
A = (Cnum2val/ Cden2val);
B = (Cnum1val/ Cden2val);
C = - (Cden1val/ Cden2val);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN : EXTRA POLE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 2 * fs * (z-1)/(z+1);
Hz = p1/(s+p1);
Hz = simplify(Hz);
[N,D] = numden(Hz);
Cnum = coeffs(N, z);
Cden = coeffs(D, z);
Cnum1 = Cnum(1);
Cnum2 = Cnum(2);
Cden1 = Cden(1);
Cden2 = Cden(2);
p1val = wp;
fsval = 2*pi*fs; % sampling frequency in the discrete time domain
Cnum1valextra = double( subs(Cnum1 , {p1 , fs},{p1val , fsval}));
Cnum2valextra = double( subs(Cnum2 , {p1 , fs},{p1val , fsval}));
Cden1valextra = double( subs(Cden1 , {p1 , fs},{p1val , fsval}));
Cden2valextra = double( subs(Cden2 , {p1 , fs},{p1val , fsval}));
Aextra = ( Cnum2valextra/ Cden2valextra);
Bextra = ( Cnum1valextra/ Cden2valextra);
Cextra = - ( Cden1valextra/ Cden2valextra);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN : INTEGRATOR %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 2 * fs * (z-1)/(z+1);
Hz = 1/s;
Hz = simplify(Hz)
[N,D] = numden(Hz);
Cnum = coeffs(N, z);
Cden = coeffs(D, z);
Cnum1 = Cnum(1);
Cnum2 = Cnum(2);
Cden1 = Cden(1);
Cden2 = Cden(2);
fsval = 2*pi*fs; % sampling frequency in the discrete time domain
Cnum1valint = double( subs(Cnum1 , fs , fsval));
Cnum2valint = double( subs(Cnum2 , fs , fsval));
Cden1valint = double( subs(Cden1 , fs , fsval));
Cden2valint = double( subs(Cden2 , fs , fsval));
Aint = ( Cnum2valint/ Cden2valint);
Bint = ( Cnum1valint/ Cden2valint);
Cint = - ( Cden1valint/ Cden2valint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Noise to be added to the input of the DCO %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foff = 1e6; % [Hz]
noise_at_foff = -130; % [dBc/Hz]
noise_var = ( foff^2)/(kv^2)*10^( noise_at_foff /10);
Ts = 1/fref;
AdditiveNoiseDCOinput = sqrt( noise_var/Ts)* randn(1,length (time_array));
VarAdditiveNoiseDCOinput = var(AdditiveNoiseDCOinput);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Dither added to the TDC input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
swichTDCdtr = 1; % set swichTDCdtr = 0 in order to switch off the TDC dither
swichDCOdtr = 1; % set swichTDCdtr = 0 in order to switch off the DCO dither
dither = swichTDCdtr *(1/3)*randn(1, length( time_array)); % TDC dither
for n = 2:1:length( time_array)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Freq to Phase conversion %%%
  %%% of the input Fref %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  phase_input(n) = (Cnum2valint/ Cden2valint) * Finput(n) + ( Cnum1valint/ Cden2valint) * Finput(n-1) - ( Cden1valint/ Cden2valint) * phase_input(n-1) ;
  if n == 2
    N = 100;
  end
  %%%%%%%%%%%%%%%
  %%% TDC %%%
  %%%%%%%%%%%%%%%
  resTDC = 10e-12;
  phase_err(n) = (phase_input (n) - output_phase_dco(n-1));
  tin(n) = (1/fref)* phase_err(n);
  tin(n) = tin(n) + resTDC*dither(n);
  outTDC(n) = resTDC*floor(tin(n)/resTDC);
  eq(n) = resTDC*(tin(n)/resTDC - floor(tin(n)/resTDC));
  outTDCwithGain(n) = kp* outTDC(n)* fref;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% FILTER : zero plus integrator %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  out(n) = Cnum2val/ Cden2val * outTDCwithGain(n) + Cnum1val/ Cden2val * outTDCwithGain(n-1) - Cden1val/ Cden2val * out(n-1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% FILTER: Extra Pole %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outExtra(n) = ( Cnum2valextra/ Cden2valextra) * out(n) + (Cnum1valextra/ Cden2valextra) * out(n-1) - (Cden1valextra/ Cden2valextra) * outExtra(n-1) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% Digital Filter Gain %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outExtraGain(n) = klf * outExtra(n);
  input_dco(n) = outExtraGain(n); % + AdditiveNoiseDCOinput(n);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% DCO with quantization %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  c = 10e-12;
  l = 100e-9;
  kdcoequivalent = -l /(2*sqrt((l*c)^3));
  deltaC = kv*(1/ kdcoequivalent) * input_dco(n);
  resDCO = 1e-20; % [Farads] 1e-20 F is negligible
  ditherDCO(n) = swichDCOdtr * floor(8* rand(1,1) -4);
  deltaCq = resDCO * (floor(deltaC/resDCO) + ditherDCO(n));
  eqDeltaC(n) = resDCO*(deltaC/resDCO - floor(deltaC/resDCO));
  freq_dco(n) = 1/ sqrt(l*c) + kdcoequivalent* deltaCq + kv*AdditiveNoiseDCOinput(n);
  output_phase_dco(n) = ( Cnum2valint/ Cden2valint) * freq_dco(n)/N + (Cnum1valint/ Cden2valint) * freq_dco(n-1)/N - (Cden1valint/ Cden2valint) *  output_phase_dco(n-1);
end
%%%%% Selection of the steady state samples %%%%%
SteadyStateOffset = 1e5;
input_dco_steady = input_dco( SteadyStateOffset:end) ;
outExtraGain_steady = outExtraGain(SteadyStateOffset: end) ;
freq_dco = freq_dco (SteadyStateOffset:end) ;
%%%%% Calculate the DCO ’s Phase Noise %%%%%
num_segments = 16;
vin = freq_dco/kv; %; input_dco;
phase = filter(Ts*2*pi*kv ,[1 -1],vin -mean(vin));
window_length = floor(length(phase)/ num_segments);
[PSDphase ,f] = pwelch(phase , window_length ,[],[] ,1/Ts ,'twosided');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculate the PSD of %%%%%
%%%%% the quantization error %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
window_length = floor(length(eq)/ num_segments);
[PSDeq ,feq] = pwelch(eq ,window_length ,[] ,[],1/Ts ,'twosided');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Analytical Preditctions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 2*pi*f/fref;
H = (A+B*exp(-1i*w))./(1-C*exp(-1i*w));
Hint = (Aint+Bint*exp(-1i*w))./(1-Cint*exp(-1i*w));
Hextra = (Aextra+Bextra*exp(-1i*w))./(1- Cextra*exp(-1i*w));
%%%% Feedforward and feedback transfer functions for the phase noise
%%%% of the reference oscillator
HffRef = kv*kp * klf * H .* Hint .* Hextra; %Hfeedforward
HfbRef = 1/N; % Hfeedback
HclRef = HffRef ./(1+ HffRef.*HfbRef); % Closed Loop TF for the phase noise of Ref
%%%% Feedforward and feedback transfer functions for noise at the DCOinput
HffDCO = Hint*kv; % Hfeedforward
HfbDCO = (1/N) * kp * klf * H .* Hextra; % Hfeedback
HclDCO = HffDCO ./(1+ HffDCO.*HfbDCO); % Closed Loop TF for the noise at the DCO input
VarTDCadditiveNoise = (fref).^2*(( swichTDCdtr* resTDC)^2) /12;
VarTDCquant = ( fref).^2*(resTDC ^2)/12;
PSDinputDCO_TDCquantComponent = (abs(HclRef)./abs(Hint*kv)).^2 * ( VarTDCquant) ;
PSDinputDCO_TDCAdditive = (abs(HclRef)./abs(Hint*kv)).^2 * (VarTDCadditiveNoise);
accumulationGain = (1/(Ts*2*pi));
% The sdomain transfer function of an accumulator that is clocked by
% a signal with period Ts is TFintegrator = ((1/Ts)/(2*pi) )/f
PSDphasePredicted_TDCquant = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*PSDinputDCO_TDCquantComponent;
PSDphasePredicted_TDCAdditive = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*PSDinputDCO_TDCAdditive;
PSDinputDCO_DCOcomponent = (abs( HclDCO./(Hint*kv))).^2 * VarAdditiveNoiseDCOinput;
PSDphasePredicted_DCOinput = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*( PSDinputDCO_DCOcomponent);
VarDCOquant = (kdcoequivalent/kv) .^2*(resDCO^2)/12;
PSDinputDCO_DCOquant = (abs( HclDCO./(Hint*kv))).^2 * VarDCOquant;
PSDphasePredicted_DCOquant = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*( PSDinputDCO_DCOquant);
PSDeqPredicted = (fref).^2*Ts*ones (1,length(PSDeq)) *((resTDC)^2)/12;
%%%%% Plot PSD Quantization error %%%%%
figure
h = semilogx( feq , 10*log10(( fref).^2*PSDeq), 'b', feq , 10*log10( PSDeqPredicted), 'r');
set(gca , 'fontsize', 15, 'fontweight', 'bold')
xlabel('Frequency [Hz]')
ylabel('PSD [dB/Hz]')
grid on
set(h(1), 'LineWidth', 6)
set(h(2), 'LineWidth', 4)
axis([1e3 0.5e7 -190 -130])
legend('Matlab Model','White Noise Approximation')
%%%%% Plot Phase Noise %%%%%
figure
h=semilogx( ...
  f, 10*log10( PSDphase), 'b',
  f, 10*log10(PSDphasePredicted_TDCAdditive +PSDphasePredicted_TDCquant +PSDphasePredicted_DCOinput + PSDphasePredicted_DCOquant ), 'r', ...
  f, 10*log10( PSDphasePredicted_DCOinput), 'g--', ...
  f, 10*log10( PSDphasePredicted_TDCquant +PSDphasePredicted_TDCAdditive), 'g');
grid on
set(gca , 'fontsize', 15, 'fontweight', 'bold')
xlabel('Frequency [Hz]')
ylabel('Phase Noise [dBc/Hz]')
grid on
set(h(1), 'LineWidth', 9)
set(h(2), 'LineWidth', 8)
set(h(3), 'LineWidth', 2)
set(h(4), 'LineWidth', 3)
axis([1e3 0.5e7 -150 -74])
legend('Matlab Model','Analytical Predictions', 'DCO noise', 'TDC noise')
