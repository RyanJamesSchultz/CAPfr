% Comparing the probability of exceedance between bound/unbound sequences of earthquakes.
% Used to make Figure S41.
clear;

% Define some constants.
Mc=0.0;
b=1.0;
SI=-2.0;
Vt=1e4;
Mmax=2.5;

% Define some dependent constants.
Nc=Vt*10^(SI-b*Mc);
M=Mc:0.01:Mmax+2.0;

% Get the catalogues.
[PDFu,CDFu,Nu]=GR_MFD(M,Mc,Inf, log10(Nc),b,'count');
[PDFb,CDFb,Nb]=GR_MFD(M,Mc,Mmax,log10(Nc),b,'count');

% Get the probabilites of exceedance.
Peu=1-exp(-Nu);
Peb=1-exp(-Nb);


% Plot.
figure(541); clf;
% Pe vs M.
subplot(211);
semilogy(M,Nu,'-m'); hold on;
semilogy(M,Nb,'-b');
semilogy([3 3],ylim,'--','Color','#691313');
semilogy([1 1],ylim,'--','Color','#EDB120');
xlabel('Magnitude (Mw)'); ylabel('Counts');
% Pe vs M.
subplot(212); hold on;
plot(M,Peu,'-m');
plot(M,Peb,'-b');
plot([3 3],ylim,'--','Color','#691313');
plot([1 1],ylim,'--','Color','#EDB120');
xlabel('Magnitude (Mw)'); ylabel('Probability of Exceedance');

