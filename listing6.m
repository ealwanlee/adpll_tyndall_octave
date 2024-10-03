close all
clear all
tstep = 1;
nsamples = 2^20;
time = 0:tstep:nsamples-1;
period = 2^9;
osc = sin(2*pi*(1/period)*time);
hannWindow = hann(round(length(time)/32));
[PSDoutweltch ,wappo] = pwelch(osc, hannWindow);
PSDoutweltch = pi*PSDoutweltch;
fnorm = 0:0.5/(length(PSDoutweltch) -1):0.5;
figure
h = semilogx(fnorm, 10*log10(PSDoutweltch), 'b','LineWidth',4);
set(gca, 'fontsize', 15, 'fontweight', 'bold')
xlabel('Normalized Frequency')
ylabel('PSD [dB/sample]')
axis([1e-5 0.5 -120 45])
grid on
