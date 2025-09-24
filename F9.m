% Script that makes an expected GR-MFD, based on a mixture of truncated/tapered GR-MFDs.
% Used to make Figure 9.
clear;

% Predefine some constants.
%BoundCases={'PNR1z-all.mat', 'St1-all.mat', 'FORGE-c1.mat','FORGE-c2a.mat','FORGE-c2.mat', 'PNR2-cW.mat'};
CaseFile='FORGE-c2a.mat';
dm2=0.45;

% Grab the data.
load(CaseFile,'D');

% Get some relevant parameters.
M=D.M(D.M>=D.Mc);
T=D.T(D.M>=D.Mc);
N=length(M);
m1=D.Mc;
b=mean(D.b);
m=m1:0.01:max(D.M+0.5);

% Fit the GR-MFD, but just to get the plotting outputs.
[~,~,a,~,~,Mgr,Ngr,ngr]=Bval(D.M,D.Mc,1.5*D.dMb);
po=[-mean(b),mean(a)];
Mgr_fit=[m1, max(M)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Get the expected Mmax values.
Vi=interp1(D.t,D.V,T,'linear','extrap');
%Vi=D.V;
m2m=Mmax_V(Vi,D.G,'McGarr');
m2g=Mmax_V(Vi,D.g,'Galis');

% Compute the expected GR-MFDs.
[PDFu,CDFu,SVFu]=GR_MFD(m,m1,Inf,a,b,'norm');
%[PDFt,CDFt,SVFt]=GR_MFDtaper(m,m1,0.7-dm2,a,b,'norm');
[PDFm,CDFm,SVFm]=GR_MFD_comp(m,m1,m2m,a,b*ones(size(m2m)));
[PDFmt,CDFmt,SVFmt]=GR_MFD_compt(m,m1,m2m-dm2,a,b*ones(size(m2m)));
[PDFg,CDFg,SVFg]=GR_MFD_comp(m,m1,m2g,a,b*ones(size(m2g)));

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#eab3fa'};
styles={'-','--','-'};
names={'McGarr','Galis','Unbound'};
GREY=[0.85,0.85,0.85];

% Plot.
figure(9); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k','HandleVisibility','off'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY,'HandleVisibility','off');
%semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'k','HandleVisibility','off');
%semilogy(m,SVFt*N,'-','LineWidth',2,'LineStyle','--','Color',colours{3},'DisplayName',['M_{MAX} taper']);
semilogy(m,SVFu*N,'-','LineWidth',2,'LineStyle',styles{3},'Color',colours{3},'DisplayName',['M_{MAX} ',names{3}]);
%semilogy(m,PDFu*10^a,'-','LineWidth',1,'LineStyle',styles{3},'Color',colours{3},'HandleVisibility','off');
semilogy(m,SVFg*N,'-','LineWidth',2,'LineStyle',styles{2},'Color',colours{2},'DisplayName',['M_{MAX} ',names{2}]);
%semilogy(m,PDFg*10^a,'-','LineWidth',1,'LineStyle',styles{2},'Color',colours{2},'HandleVisibility','off');
semilogy(m,SVFm*N,'-','LineWidth',2,'LineStyle',styles{1},'Color',colours{1},'DisplayName',['M_{MAX} ',names{1}]);
semilogy(m,SVFmt*N,'-','LineWidth',2,'LineStyle','--','Color',colours{1},'DisplayName',['M_{MAX} ',names{1},' taper']);
%semilogy(m,PDFm*10^a,'-','LineWidth',1,'LineStyle',styles{1},'Color',colours{1},'HandleVisibility','off');
xlim([min(Mgr)-D.dMd/2 max(m)+D.dMd/2]); ylim([0.7*1e0 1.3*max(Ngr)]);
plot(m1*[1 1],ylim,'--k','HandleVisibility','off');
xlabel('Magnitude'); ylabel('Count');
legend('Location','northeast');

