%%%%%%%% Jitter extraction %%%%%%%%%%%
%%%%%% Processing the Periods %%%%%%%
idxCrossingsPositive = [];
idxCrossingsNegative = [];
for idx = 2:1:nsamples
  %%%% rising edge detection %%%%
  if (sinsignalNoisy(idx) > 0 ) && ( sinsignalNoisy(idx-1) < 0 )
    idxCrossingsPositive = [idxCrossingsPositive idx];
    idxCrossingsNegative = [idxCrossingsNegative (idx-1)];
  end
end
grid on
%%%% Interpolation %%%%%
xa = time(idxCrossingsNegative);
xb = time(idxCrossingsPositive);
ya = sinsignalNoisy(idxCrossingsNegative);
yb = sinsignalNoisy(idxCrossingsPositive);
te = xb - (xb - xa)/(yb - ya)*yb; % crossings instants
teideal = period*(1:1:length(te)); % ideal crossings
deltaT = (teideal - te);
deltaPhi = 2*pi*deltaT./(period); % Calculate DeltaPhi from DeltaT
jitterExtracted = deltaPhi;
figure
h=plot(
  1:1:length(jitter), jitter, '-bo', ...
  (1:1:length(jitterExtracted))/length(jitterExtracted)*length(jitter),(jitterExtracted), 'r.' ...
);
set(h(1), 'LineWidth', 4);
set(h(2), 'MarkerSize', 5);
set(gca, 'fontsize', 15, 'fontweight', 'bold')
xlabel('Time Index')
ylabel('Phase Deviation')
legend('Actual','Extracted')
grid on
