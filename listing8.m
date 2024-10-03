hannWindow = hann(round(length(jitterExtracted)/32));
[PSDjitterExtracted ,wappo] = pwelch(jitterExtracted, hannWindow);
PSDjitterExtracted = pi*PSDjitterExtracted;
Nperiods = nsamples/period;
SamplesPerPeriod = nsamples/Nperiods;
fjitterExtracted =(0:(0.5/SamplesPerPeriod)/(length(PSDjitterExtracted)-1):(0.5/SamplesPerPeriod)) ;
fjitterExtracted = fjitterExtracted + 1/period;
figure
h = semilogx(fJitter, 10*log10(PSDphaseComponent/2),'b', fjitterExtracted, 10*log10(SamplesPerPeriod*PSDjitterExtracted/2), 'g');
set(h(1), 'LineWidth', 6);
set(h(2), 'LineWidth', 2);
set(gca, 'fontsize', 15, 'fontweight', 'bold')
xlabel('Normalized Frequency')
ylabel('PSD [dB/sample]')
legend('Phase Deviation','Extracted Phase Deviation')
axis([(1/period)-3e-4 0.5 -120 45])
grid on
