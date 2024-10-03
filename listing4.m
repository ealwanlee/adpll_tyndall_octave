clc
close all
clear all
tstep = 1;
nsamples = 2e5; %2^23;
time = 0:tstep:nsamples-1;
period = 2^9;
ThreeSigma = period/20;
Sigma = ThreeSigma/3;
varjitter = Sigma^2;
Nperiods = length(time)/period +10;
jitterWhite = Sigma*randn(1,nsamples);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILTER DESIGN: Single POLE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z z1 outz inz p1;
fs = tstep/period;
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
fsval = 2*pi*fs; % sampling frequency in the discrete time domain
p1val = fsval/100; % frequency of the pole
Cnum1valextra = double(subs(Cnum1, {p1, fs},{p1val,fsval}));
Cnum2valextra = double(subs(Cnum2, {p1, fs},{p1val,fsval}));
Cden1valextra = double(subs(Cden1, {p1, fs},{p1val,fsval}));
Cden2valextra = double(subs(Cden2, {p1, fs},{p1val,fsval}));
Aextra = (Cnum2valextra/ Cden2valextra);
Bextra = (Cnum1valextra/ Cden2valextra);
Cextra = - (Cden1valextra/ Cden2valextra);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FILTER: Implementation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jitterFiltered = zeros(1,length(jitterWhite));
for n = 2:1:length(jitterWhite)
  jitterFiltered(n) = (Cnum2valextra/ Cden2valextra) *jitterWhite(n) + (Cnum1valextra/ Cden2valextra) *jitterWhite(n-1) - (Cden1valextra/ Cden2valextra) *jitterFiltered(n-1) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SECOND FILTER: Design and Implementation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jitterWhiteNew = Sigma*randn(1,nsamples);
syms z z1 outz inz p1;
fs = tstep/period;
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
fsval = 2*pi*fs; % sampling frequency in the discrete time domain
p1val = fsval/10000; % frequency of the pole (Hz)
Cnum1valextra = double(subs(Cnum1, {p1, fs},{p1val,fsval}));
Cnum2valextra = double(subs(Cnum2, {p1, fs},{p1val,fsval}));
Cden1valextra = double(subs(Cden1, {p1, fs},{p1val,fsval}));
Cden2valextra = double(subs(Cden2, {p1, fs},{p1val,fsval}));
Aextra = (Cnum2valextra/ Cden2valextra);
Bextra = (Cnum1valextra/ Cden2valextra);
Cextra = - (Cden1valextra/ Cden2valextra);
jitterFilteredNew = zeros(1,length(jitterWhite));
for n = 2:1:length(jitterWhite)
  jitterFilteredNew(n) = (Cnum2valextra/ Cden2valextra)  * jitterWhiteNew(n) + (Cnum1valextra/  Cden2valextra) * jitterWhiteNew(n-1) -  (Cden1valextra/ Cden2valextra) *  jitterFilteredNew(n-1) ;
end
jitter = 0.01*( 0.01*cumsum(jitterFilteredNew) + 0.01*jitterFiltered );
sinsignalNoisy = sin(2*pi*(1./period).*time+ jitter);
