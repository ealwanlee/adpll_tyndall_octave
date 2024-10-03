hannWindow = hann(round(length(time)/16));
[PSDoutweltch ,wappo] = pwelch(sinsignalNoisy,hannWindow);
PSDoutweltch = pi*PSDoutweltch;
f = 0:0.5/(length(PSDoutweltch) -1):0.5;
hannWindow = hann(round(length(jitter)/16));
[PSDphaseComponent ,wappo] = pwelch(jitter,hannWindow);
PSDphaseComponent = pi*PSDphaseComponent;
fJitter = (0:0.5/(length(PSDphaseComponent) -1):0.5) ;
fJitter = fJitter + 1/period;
fMirror = 0:0.5/(length(PSDoutweltch) -1):1/period;
temp = PSDphaseComponent(1:1:length(fMirror));
PSDphaseComponentMirror = temp(end:-1:1);
figure
h = semilogx(
  f, 10*log10(PSDoutweltch), 'b', ...
  fJitter,10*log10(PSDphaseComponent/2), 'r', ...
  fMirror, 10*log10(PSDphaseComponentMirror/2), 'r' );
set(h(1), 'LineWidth', 6);
set(h(2), 'LineWidth', 2);
set(h(3), 'LineWidth', 2);
set(gca, 'fontsize', 15, 'fontweight', 'bold')
xlabel('Normalized Frequency')
ylabel('PSD [dB/sample]')
legend('Amplitude','Phase Deviations')
axis([1e-5 0.5 -120 45])
grid on
