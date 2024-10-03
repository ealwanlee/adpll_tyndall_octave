%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PLL MODEL WHICH INCLUDES %%%%
%%% A TDC AND A DCO WITH %%%%
%%% MACHINE PRECISION %%%%
%%% RESOLUTIONS %%%%
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
%Nsamples = 2^10;
Nsamples = 2e5;
fz = 2.2e4;
wz = 2*pi*fz;
fp = 1.734e5;
wp =2*pi*fp;
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
freq_freerunning = fref*N; % freerunning frequency of the controlled oscillator
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
eq = zeros(Nsamples ,1);
tin = zeros(Nsamples ,1);
outTDC = zeros( Nsamples ,1);
outTDCwithGain = zeros(Nsamples ,1);
eqDeltaC = zeros( Nsamples ,1);
Finput = fref * ones(Nsamples ,1); % Input Frequency
time_array =(1/fref) * (0:1: Nsamples -1); % Time Array
fs = fref; % Sampling Frequency of the Digital Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN : ZERO plus INTEGRATOR %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z z1 outz inz p1;
s = 2 * fs * (z-1)/(z+1);
Hz = (s + z1) / (s);
%Hz = simple(Hz);
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
Cnum1val = double( subs(Cnum1 , {z1 , fs},{z1val ,fsval}));
Cnum2val = double( subs(Cnum2 , {z1 , fs},{z1val ,fsval}));
Cden1val = double( subs(Cden1 , {z1 , fs},{z1val ,fsval}));
Cden2val = double( subs(Cden2 , {z1 , fs},{z1val ,fsval}));
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
Cnum1valextra = double( subs(Cnum1 , {p1 , fs},{p1val ,fsval}));
Cnum2valextra = double( subs(Cnum2 , {p1 , fs},{p1val ,fsval}));
Cden1valextra = double( subs(Cden1 , {p1 , fs},{p1val ,fsval}));
Cden2valextra = double( subs(Cden2 , {p1 , fs},{p1val ,fsval}));
Aextra = ( Cnum2valextra/ Cden2valextra);
Bextra = ( Cnum1valextra/ Cden2valextra);
Cextra = - ( Cden1valextra/ Cden2valextra);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN : INTEGRATOR %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 2 * fs * (z-1)/(z+1);
Hz = 1/s;
Hz = simplify(Hz);
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
	phase_err(n) = (phase_input (n) - output_phase_dco(n-1));
	tin(n) = (1/fref)* phase_err(n);
	resTDC = 1e-19; % [seconds] rest = 1e-19 introduces negligible quant error
	outTDC(n) = resTDC*floor(tin(n)/resTDC);
	eq(n) = resTDC* (tin(n)-floor(tin(n)));
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
	input_dco(n) = outExtraGain(n); % +AdditiveNoiseDCOinput(n);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% DCO with quantization %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	c = 10e-12;
	l = 100e-9;
	kdcoequivalent = -l /(2*sqrt((l*c)^3));
	deltaC = kv*(1/ kdcoequivalent) * input_dco(n);
	resDCO = 1e-20; % [Farads] 1e-20 F is negligible
	deltaCq = resDCO * floor(deltaC/resDCO);
	eqDeltaC(n) = resDCO*(deltaC/resDCO - 	floor(deltaC/resDCO));
	freq_dco(n) = 1/ sqrt(l*c) + kdcoequivalent* deltaCq + kv*AdditiveNoiseDCOinput(n);
	output_phase_dco(n) = ( Cnum2valint/ Cden2valint) * freq_dco(n)/N + (Cnum1valint/ Cden2valint) * 	freq_dco(n-1)/N - (Cden1valint/ Cden2valint) * output_phase_dco(n-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Selection of the steady state samples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SteadyStateOffset = 1e5;
freq_dco_steady = freq_dco( SteadyStateOffset:end) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculate the DCO ’s Phase Noise %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_segments = 16;
DCOin = freq_dco_steady/kv; %; input_dco;
phase = filter(Ts*2*pi*kv ,[1 -1], DCOin - mean(DCOin));
window_length = floor(length(phase)/ num_segments);
size(phase)
[PSDphase ,f] = pwelch(phase , window_length ,[],[] ,1/Ts ,'twosided');
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
VarTDCquant = (pi^2)*( fref).^2*(resTDC ^2)/12;
PSDinputDCO_TDCquantComponent = (abs(HclRef./(Hint*kv))).^2 * VarTDCquant;
accumulationGain = (1/(Ts*2*pi));
% The sdomain transfer function of an accumulator that is clocked by
% a signal with period Ts is TFintegrator = ((1/Ts)/(2*pi) )/f
PSDphasePredicted_TDCquant = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*PSDinputDCO_TDCquantComponent;
PSDinputDCO_DCOcomponent = (abs( HclDCO./(Hint*kv))).^2 * VarAdditiveNoiseDCOinput;
PSDphasePredicted_DCOinput = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*( PSDinputDCO_DCOcomponent);
VarDCOquant = (kdcoequivalent/kv) .^2*(resDCO^2)/12;
PSDinputDCO_DCOquant = (abs( HclDCO./(Hint*kv))).^2 * VarDCOquant;
PSDphasePredicted_DCOquant = (1./f).^2.* (Ts*2*pi*kv)^2 *accumulationGain ^2* Ts.*( PSDinputDCO_DCOquant);
%%%%% Plot Phase Noise %%%%%
figure
h=semilogx( f, 10*log10( PSDphase), 'b', f, 10*log10(PSDphasePredicted_TDCquant + PSDphasePredicted_DCOinput + PSDphasePredicted_DCOquant ), 'r');
grid on
set(gca , 'fontsize', 15, 'fontweight','bold')
xlabel('Frequency [Hz]')
ylabel('Phase Noise [dBc/Hz]')
grid on
set(h(1), 'LineWidth', 2)
set(h(2), 'LineWidth', 2)
axis([1e4 0.5e7 -150 -104])
legend('Matlab Model','Analytical Predictions')
