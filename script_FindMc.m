% Script to optimize/select the magnitude-of-completeness for a case.
clear;

% Define some constants.
dMc=-0.50:0.001:+0.50;
%CaseList={'FORGE22','FORGE24', 'PNR2','PNR1z', 'St1', 'SSFS93','SSFS00','SSFS03','SSFS04','SSFS05'};
CaseList={'SSFS93'};
i=[1:2];

% Get all of the case data.
D=PreProcData(CaseList,i);

% Loop over all of the case data.
for k=1:length(D)
    
    % Get the case's catalogue.
    Lat=D(k).Lat'; Lon=D(k).Lon'; Dep=D(k).Dep';
    T=D(k).T'; M=D(k).M';
    m1b=D(k).Mc;
    m1k=D(k).Mk;
    dMb=D(k).dMb;
    dMd=D(k).dMd;
    
    % Loop over the bootstrap trials.
    for i=1:length(dMc)
        i
        
        % Apply a perturbation.
        Mi=M;
        m1bi=m1b+dMc(i);
        
        % Truncate on the magnitude of completeness.
        Imb=(Mi>=m1bi);
        
        % Derive some variables.
        Ncb=length(Mi(Imb));
        
        % Fit the GR-MFD b-value.
        [bi,b_err,ai,R2i,LL,Mgr,Ngr,ngr]=Bval(Mi,m1bi,dMb);
        b(i)=bi;
        a(i)=ai;
        R2(i)=R2i;
        nLL(i)=-GR_MFD_LL(Mi(Imb),m1bi,Inf,ai,bi);
        
        % Get the expected Mlrg value.
        Mlrg_est(i)=m1bi+log10(Ncb)/bi;
        
    end
    
    % Save data into the output structure.
    D(k).b=b;
    D(k).R2=R2;
    D(k).nLL=nLL;
    D(k).Mlrg_est=Mlrg_est;
    
end


% Plot.
figure(4); clf;
nm=round(1.5*D(k).dMd/mean(diff(dMc)));
% b-values.
subplot(311);
%plot(dMc,D(k).b); hold on;
plot(dMc,movmean(D(k).b,nm));
xlabel('Mc Perturbation'); ylabel('b-value');
% R2 values.
subplot(312);
%plot(dMc,D(k).R2); hold on;
plot(dMc,movmean(D(k).R2,nm));
xlabel('Mc Perturbation'); ylabel('R^2 value');
% nLL values.
subplot(313);
%plot(dMc,D(k).nLL); hold on;
plot(dMc,movmean(D(k).nLL,nm));
xlabel('Mc Perturbation'); ylabel('nLL value');
