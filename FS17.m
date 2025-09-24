% GR-MFD for the clustered catalogues at FORGE.
% Used to make Figure S17.
clear;

% Define some constants.
dMb=0.1;
dMd=0.1;

% Get all of the case data.
D1 =PreProcData({'FORGE22'},1:2); % Cluster 1.
D2a=PreProcData({'FORGE22'},  3); % Cluster 2a.
D2 =PreProcData({'FORGE24'},1:4); % Cluster 2.
D3 =PreProcData({'FORGE24'},5:8); % Cluster 3.
D(1)=D1; D(2)=D2a; D(3)=D2; D(4)=D3; 

% Loop over all of the k case data.
for k=1:length(D)
    
    % Get the case's catalogue.
    Lat=D(k).Lat'; Lon=D(k).Lon'; Dep=D(k).Dep';
    T=D(k).T'; M=D(k).M';
    m1b=D(k).Mc;
    m1k=D(k).Mk;
    
    % Fit the GR-MFD.
    [b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(M,m1b,dMb);
    
    % Get the GR-MFD plotting values.
    po=[-mean(b),a];
    Mgr_fit=[m1b, max(M)];
    Ngr_fit=10.^polyval(po,Mgr_fit);
    
    % Get the expected Mlrg value.
    Ncb=length(M(M>=m1b));
    Mlrg_est=m1b+log10(Ncb)/b;
    
    % Report b-values and observed-expected Mlrg discrepancy.
    [b,b_err]
    max(M)-Mlrg_est
    R2
    
    % Save data into the output structure.
    D(k).b=b;
    D(k).Mgr=Mgr;
    D(k).Ngr=Ngr;
    D(k).ngr=ngr;
    D(k).Mgr_fit=Mgr_fit;
    D(k).Ngr_fit=Ngr_fit;
    D(k).Mlrg_est=Mlrg_est;
    
end


% Plot catalogue filtering info.
GREY=[0.85,0.85,0.85];
figure(517); clf;
% GR-FMD C1.
subplot(221);
semilogy(D(1).Mgr, D(1).Ngr, 'o', 'Color', 'k'); hold on;
bar(D(1).Mgr,D(1).ngr, 'FaceColor', GREY);
semilogy(D(1).Mgr_fit, D(1).Ngr_fit, '-', 'Color', 'black');
xlim([min(D(1).Mgr)-dMd/2 max(D(1).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(1).Ngr)]);
plot(D(1).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (Mw)'); ylabel('Count');
% GR-FMD C2a.
subplot(222);
semilogy(D(2).Mgr, D(2).Ngr, 'o', 'Color', 'k'); hold on;
bar(D(2).Mgr,D(2).ngr, 'FaceColor', GREY);
semilogy(D(2).Mgr_fit, D(2).Ngr_fit, '-', 'Color', 'black');
xlim([min(D(2).Mgr)-dMd/2 max(D(2).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(2).Ngr)]);
plot(D(2).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (Mw)'); ylabel('Count');
% GR-FMD C2.
subplot(223);
semilogy(D(3).Mgr, D(3).Ngr, 'o', 'Color', 'k'); hold on;
bar(D(3).Mgr,D(3).ngr, 'FaceColor', GREY);
semilogy(D(3).Mgr_fit, D(3).Ngr_fit, '-', 'Color', 'black');
xlim([min(D(3).Mgr)-dMd/2 max(D(3).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(3).Ngr)]);
plot(D(3).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (Mw)'); ylabel('Count');
% GR-FMD C3.
subplot(224);
semilogy(D(4).Mgr, D(4).Ngr, 'o', 'Color', 'k'); hold on;
bar(D(4).Mgr,D(4).ngr, 'FaceColor', GREY);
semilogy(D(4).Mgr_fit, D(4).Ngr_fit, '-', 'Color', 'black');
xlim([min(D(4).Mgr)-dMd/2 max(D(4).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(4).Ngr)]);
plot(D(4).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (Mw)'); ylabel('Count');
% GR-FMD C4.
%subplot(325);
%semilogy(D(5).Mgr, D(5).Ngr, 'o', 'Color', 'k'); hold on;
%bar(D(5).Mgr,D(5).ngr, 'FaceColor', GREY);
%semilogy(D(5).Mgr_fit, D(5).Ngr_fit, '-', 'Color', 'black');
%xlim([min(D(5).Mgr)-dMd/2 max(D(5).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(5).Ngr)]);
%plot(D(5).Mc*[1 1],ylim,'--k');
%xlabel('Magnitude (Mw)'); ylabel('Count');

