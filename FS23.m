% GR-MFD for the full catalogue at PNR-2.
% Used to make Figure S23.
clear;

% Define some constants.
dMb=0.1;
dMd=0.1;

% Get all of the case data.
D=PreProcData({'PNR2'},1:7);

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
figure(523); clf;
% PNR-2 2019.
semilogy(D(1).Mgr, D(1).Ngr, 'o', 'Color', 'k'); hold on;
bar(D(1).Mgr,D(1).ngr, 'FaceColor', GREY);
semilogy(D(1).Mgr_fit, D(1).Ngr_fit, '-', 'Color', 'black');
xlim([min(D(1).Mgr)-dMd/2 max(D(1).Mgr)+dMd/2]); ylim([0.7 1.3*max(D(1).Ngr)]);
plot(D(1).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (Mw)'); ylabel('Count');


