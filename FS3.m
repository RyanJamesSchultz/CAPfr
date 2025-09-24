% Script to perform CAP-tests on simulated data.
% Used to make Figures S3-S6.
clear;

% Define some simulation constants.
Nc=1e3;
m1=0.0;
b=1.0;
d=0.01;

% Define some CAP-test constants.
Ns=1e1;
Rflag1='resample';
Rflag2='none';
Nr1=1e2;
Nr2=0;
SMOOTHflag='mean';
Kgr=2+1;
dMm=0.1;
q=-1;


% Define a hypothetical injection scenario.
t=1:2000; % mins.
v=zeros(size(t)); v(10:200)=linspace(0,10,191); v(201:1500)=10; % m3/min.
V=cumsum(v);

% Get the timing of events.
Vi=linspace(max(V)/Nc,max(V),Nc);
T=interp1(V(v>0),t(v>0),Vi,'linear');

% Make the true Mmax.
%m2=linspace(Inf,Inf,Nc); scenario='Unbound';  % FS3.
%m2=linspace(2.3,2.3,Nc); scenario='Tectonic'; % FS4.
%m2=Mmax_V(Vi,1e-1,'McGarr'); scenario='McGarr'; % FS5.
m2=Mmax_V(Vi,5e+6, 'Galis'); scenario='Galis'; % FS6.
[m2(1) m2(end)]

% Get a catalogue.
a=log10(Nc);
M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);

% Loop over all of the case data.
for k=1:1
    
    % Get the case's catalogue.
    %Lat=D(k).Lat'; Lon=D(k).Lon'; Dep=D(k).Dep';
    %T=D(k).T'; M=D(k).M';
    m1b=m1;
    m1k=m1;
    dMb=0.1;
    dMd=0.1;
    
    % Loop over the bootstrap trials.
    for i=1:Ns
        i
        
        % Apply a perturbation.
        Mi=M+(dMd*rand(size(M))-dMd/2);
        m1bi=m1b+((dMd/2)*rand(size(m1b))-dMd/4);
        m1ki=m1k+((dMd/2)*rand(size(m1k))-dMd/4);
        
        % Truncate on the magnitude of completeness.
        Imb=(Mi>=m1bi);
        Imk=(Mi>=m1ki);
        
        % Derive some variables.
        Ncb=length(Mi(Imb));
        Nck=length(Mi(Imk));
        
        % Fit the GR-MFD b-value.
        [bi,b_err,ai,R2,~,Mgr,Ngr,ngr]=Bval(Mi,m1bi,dMb);
        b(i)=bi;
        a(i)=ai;
        
        % Get the expected Mlrg value.
        Mlrg_est(i)=m1bi+log10(Ncb)/bi;
        
        % Do the MLE-fitting.
        [m2_fit,llr]=M2fit(Mi(Imb),m1bi,bi,Nr1);
        m2_avg(i)=mean(m2_fit);
        m2_std(i)=std(m2_fit);
        m2_p90(i)=prctile(m2_fit,90);
        
        % Do the KS-test (and regular KS-test).
        Msamp=GR_MFD_Rand(m1bi,m2_avg(i), log10(Nck),bi, [1 Nck]);
        [~,p_full]=kstest2(Mi(Imb),Msamp);
        KSp_fl(i)=p_full;
        KSp_dm(i)=KS_dM_test(Mi(Imk),m1ki,Rflag1);
        
    end
    
    % Save data into the output structure.
    D(k).b=b;
    D(k).Mlrg_est=Mlrg_est;
    D(k).m2=m2_fit;
    D(k).KSp_fl=KSp_fl;
    D(k).KSp_dm=KSp_dm;
    
    % Truncate on the magnitude of completeness.
    Imb=(M>=m1b);
    Ncb=length(M(Imb));
    
    % Get the sequence of largest events (and indicies to when they occur).
    Mlrg=OrderStatistic(M(Imb),Ncb-0,'none');
    Tlrg=T(Imb);
    
    % Interpolate
    [t,I]=unique(t);
    V=V(I); v=v(I);
    Vi=interp1(t,V,Tlrg,'linear','extrap');
    
    % Fit the McGarr Mmax model.
    Afxn = @(G) min(Mmax_V(Vi,G,'McGarr')-Mlrg)-dMm;
    G=fzero(Afxn,[1e-10 1e10]);
    Mmax_M=Mmax_V(Vi,G,'McGarr');
    
    % Fit the Galis Mmax model.
    Afxn = @(g) min(Mmax_V(Vi,g,'Galis')-Mlrg)-dMm;
    g=fzero(Afxn,[1e-10 1e15]);
    Mmax_G=Mmax_V(Vi,g,'Galis');
    
    % Regular Mmax models.
    %Mmax_M=Mmax_V(Vi,30,'McGarr');
    %Mmax_G=Mmax_V(Vi,1.5e9,'Galis');
    
    % Propose possible three possible models of Mmax.
    Sm(1).Mmax=Mmax_M;
    Sm(1).K=Kgr+1;
    Sm(2).Mmax=Mmax_G;
    Sm(2).K=Kgr+1;
    Sm(3).Mmax=(max(M)+dMm)*ones(size(Mlrg));
    Sm(3).K=Kgr+1;
    Sm(4).Mmax=(Inf)*ones(size(Mlrg));
    Sm(4).K=Kgr+0;
    
    % Do the EW-test.
    W=EnsembleW(M(Imb),m1b,Sm,mean(b),Nr2,Rflag2,SMOOTHflag);
    
    % Get each Mmax model's estimate of the NLE's magnitude.
    MnleE=zeros(size(Mlrg)); Wb=[];
    for j=1:length(W)
        W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,mean(b),q);
        MnleE=MnleE+(W(j).Mnle.*W(j).W(2:end));
        Wb=[Wb;W(j).W];
    end
    
    % Get the odds ratios, relative to the unbound model.
    OR=Wb./Wb(end,:);

    % Save some more data into the output structure.
    D(k).G=G;
    D(k).g=g;
    D(k).W=W;
    D(k).Wb=Wb;
    
    % Get the GR-MFD plotting values.
    po=[-mean(b),mean(a)];
    Mgr_fit=[m1b, max(M)];
    Ngr_fit=10.^polyval(po,Mgr_fit);
    
    % Report the case name
    disp('CASE');
    scenario
    
    % Report the simple-test results.
    disp('Simple-tests');
    [max(M) mean(Mlrg_est) max(M)-mean(Mlrg_est)]
    (1-10.^((mean(Mlrg_est)-max(M))*mean(b))/Ncb)^Ncb
    Ncb
    
    % Report the KS-test p-value.
    disp('KS(ΔM)-tests');
    geomean(KSp_dm)
    geomean(KSp_fl)
    
    % Report the MLE-fitted Mmax (and related) values.
    disp('MLE(ΔM)-tests');
    [mean(m2_avg) mean(m2_std)]
    
    % Report the model weights and the relative odds of each model.
    disp('EW(ΔM)-tests');
    Wb(:,end)
    [Wb(:,end)/Wb(end,end)]
    
end

% Save data.
save('CaseData_temp.mat','D');





% Plot.
%%
% Define some colours I'd like to use.
colours={'#345da7','#587aff','#77aaff','#eab3fa'};
styles={'-','--',':','-'};
names={'McGarr','Galis','Tectonic','Unbound'};
GREY=[0.85,0.85,0.85];
Rs=getMscale(M)*200; Rs=Rs*(100/max(Rs));

% Define expected Mlrg values.
n=1:Ncb;
Ml=log10(1:Ncb)/mean(b)+m1b;
ql=0.10; qm=0.50; qh=0.90;
Ml1=Ml-log10(n.*(1-ql.^(1./n)))/mean(b);
Ml2=Ml-log10(n.*(1-qm.^(1./n)))/mean(b);
Ml3=Ml-log10(n.*(1-qh.^(1./n)))/mean(b);

% Plot time series data.
figure(1); clf;
f=0.95*(max(M)-m1)/max(v);
% MvT plot.
scatter(T,M,Rs,'r','filled','HandleVisibility','off'); hold on;
plot(T,cummax(M),'-r');
plot(t,v*f+m1,'-b');
plot(xlim(),m1*[1 1],'--r');
plot(t,v*f+m1,'-b');
xlabel('Time'); ylabel('Magnitude');
ylim([m1-0.2 max(M)]);

% Plot EW-test results.
figure(2); clf;
% MvT plot.
ax1=subplot(311);
scatter([0,n]+1,[m1b,M(Imb)],[min(Rs(Imb)),Rs(Imb)],'r','filled','HandleVisibility','off'); hold on;
plot([0,n]+1,[m1b,Ml2] ,'-.r','HandleVisibility','off');
plot([0,n]+1,[m1b,Ml1] ,':r','HandleVisibility','off');
plot([0,n]+1,[m1b,Ml3] ,':r','HandleVisibility','off');
plot([0,n]+1,[m1b,Mlrg],'-r','HandleVisibility','off');
for i=1:length(W)
    plot([0,n]+1,[W(i).Mmax(1),W(i).Mmax],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{MAX} ',names{i}]);
    %plot([0,n]+1,[W(i).Mnle(1),W(i).Mnle],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{NLE} ',names{i}]);
end
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 Ncb+1]); ylim([m1b,max([max(Ml2),max(Mlrg)]+0.3)]);
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvT plot.
ax2=subplot(312);
bp=bar([0,n]+1.5,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(W)
    plot([1 Ncb+1],[i i]/length(W),'-w');
    bp(i).FaceColor=colours{i};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1,Ncb+1]); ylim([0 1]);
set(gca, 'XScale', 'log');
% ORvN plot.
ax3=subplot(313);
for i=1:length(W)-1
    loglog([0,n]+1,OR(i,:),'-','Color',colours{i}, 'DisplayName',['M_{MAX} ',names{i}]); hold on;
end
loglog(xlim(),1*[1 1],'-k');
loglog(xlim(),3*[1 1],':k');
loglog(xlim(),10*[1 1],':k');
loglog(xlim(),100*[1 1],':k');
xlabel('Chronological Event Number'); ylabel('Relative Odds Ratio');
xlim([1 Ncb+1]); ylim([1e-1 1e+3]);
linkaxes([ax1 ax2 ax3],'x');






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end