% Show locations of earthquakes and EW-test results for PNR1z.
% Used to make Figure 3.
clear;

% Predefine some values.
Clims=[-2 +2];

% Define an input structure.
S=struct('Stage',[],'Ploc',[]);

% Define the perf/stage locations.
S( 1).Stage=1;  S( 1).Ploc(1,:)=[3374 3376];  % S1.
S( 2).Stage=2;  S( 2).Ploc(1,:)=[3356 3358];  % S2.
S( 3).Stage=3;  S( 3).Ploc(1,:)=[3338 3340];  % S3.
S( 4).Stage=12; S( 4).Ploc(1,:)=[3180 3182];  % S12.
S( 5).Stage=13; S( 5).Ploc(1,:)=[3162 3164];  % S13.
S( 6).Stage=14; S( 6).Ploc(1,:)=[3144 3146];  % S14.
S( 7).Stage=18; S( 7).Ploc(1,:)=[3074 3076];  % S18.
S( 8).Stage=22; S( 8).Ploc(1,:)=[2997 3000];  % S22.
S( 9).Stage=30; S( 9).Ploc(1,:)=[2857 2859];  % S30.
S(10).Stage=31; S(10).Ploc(1,:)=[2840 2842];  % S31.
S(11).Stage=32; S(10).Ploc(1,:)=[2822 2824];  % S32.
S(12).Stage=37; S(10).Ploc(1,:)=[2734 2736];  % S37.
S(13).Stage=38; S(10).Ploc(1,:)=[2716 2718];  % S38.
S(14).Stage=39; S(10).Ploc(1,:)=[2699 2701];  % S39.
S(15).Stage=40; S(10).Ploc(1,:)=[2681 2683];  % S40.
S(16).Stage=41; S(10).Ploc(1,:)=[2663 2665];  % S41.

% Get all of the case data.
load('PNR1z-all.mat','D'); %D1=D;   % Cluster 1.

% Get some additional information.
for i=1:length(D)
    
    % Truncate on the magnitude of completeness.
    Imb=(D(i).M>=D(i).Mc);
    D(i).Imb=Imb;
    Ncb=sum(Imb);
    
    % Get the sequence of largest events (and indicies to when they occur).
    D(i).Mlrg=OrderStatistic(D(i).M(Imb),Ncb-0,'none');
    D(i).Tlrg=D(i).T(Imb);
    
    % Get the (log10) odds ratios, relative to the unbound model.
    OR=log10(D(i).Wb./D(i).Wb(end,:));
    D(i).OR=max(OR(1:end-1,:));
    
    % Define expected Mlrg values.
    n=1:Ncb;
    Ml=log10(1:Ncb)/mean(D(i).b)+D(i).Mc;
    ql=0.10; qm=0.50; qh=0.90;
    Ml1=Ml-log10(n.*(1-ql.^(1./n)))/mean(D(i).b);
    Ml2=Ml-log10(n.*(1-qm.^(1./n)))/mean(D(i).b);
    Ml3=Ml-log10(n.*(1-qh.^(1./n)))/mean(D(i).b);

    % Save these values.
    D(i).Ncb=Ncb;
    D(i).n=n;
    D(i).Ml1=Ml1;
    D(i).Ml2=Ml2;
    D(i).Ml3=Ml3;
    
    % Event sizes for plotting.
    Mrs=getMscale(max(vertcat(D.M)));
    D(i).Rs=getMscale(D(i).M);
    D(i).Rs=D(i).Rs*(100/Mrs);
    
    % Set the well surface location as the origin.
    D(i).Lon=D(i).Lon-D(i).Wlon(1);
    D(i).Lat=D(i).Lat-D(i).Wlat(1);
    D(i).Wlon=D(i).Wlon-D(i).Wlon(1);
    D(i).Wlat=D(i).Wlat-D(i).Wlat(1);
end

% Report the weights and odds ratio, as part of the text.
n=100+1;
w=unique([D(end).Wb]','rows','stable');
w(n,:)
max(w(n,1:end-1))/w(n,end)
log10(max(w(n,1:end-1))/w(n,end))




% Plot.

% Define some colours I'd like to use.
colours={'#345da7','#587aff','#eab3fa'};
styles={'-','--','-'};
names={'McGarr','Galis','Unbound'};
GREY=[0.85,0.85,0.85];

% Figure 3.
figure(3); clf;
% Map.
subplot(2,4,[1 2 5 6 ]); hold on;
for i=1:length(D)
    scatter(D(i).Lon(D(i).Imb), D(i).Lat(D(i).Imb), D(i).Rs(D(i).Imb), D(i).OR(2:end),'filled','MarkerEdgeColor','k');
end
plot(D(1).Wlon,   D(1).Wlat,'-k');
plot(D(1).Wlon(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    x=interp1(D(1).Wmd,D(1).Wlon,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(x,y,'-r','LineWidth',1);
    x=interp1(D(1).Wmd,D(1).Wlon,mean(S(i).Ploc,2),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,mean(S(i).Ploc,2),'linear');
    plot(x,y,'xr','LineWidth',1);
end
h = colorbar(); colormap(gca,R_colormap('Odds Ratio')); clim(Clims);
xlabel('Easting (m)'); ylabel('Northing (m)');
ylabel(h, 'Log_{10} Odds Ratio');
axis equal;
ylim([-400 200]); xlim([-1900 -800]);
% MvT.
j=1;
ax1=subplot(2,4,3:4); hold on;
scatter([0,D(j).n]+1,[D(j).Mc;D(j).M(D(j).Imb)],[min(D(j).Rs(D(j).Imb));D(j).Rs(D(j).Imb)],D(j).OR,'filled','HandleVisibility','off'); hold on;
plot([0,D(j).n]+1,[D(j).Mc,D(j).Ml2] ,'-.r','HandleVisibility','off');
plot([0,D(j).n]+1,[D(j).Mc,D(j).Ml1] ,':r','HandleVisibility','off');
plot([0,D(j).n]+1,[D(j).Mc,D(j).Ml3] ,':r','HandleVisibility','off');
plot([0,D(j).n]+1,[D(j).Mc;D(j).Mlrg],'-r','HandleVisibility','off');
for i=1:length(D(j).W)
    plot([0,D(j).n]+1,[D(j).W(i).Mmax(1),D(j).W(i).Mmax],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{MAX} ',names{i}]);
    %plot([0,D(j).n]+1,[D(j).W(i).Mnle(1),D(j).W(i).Mnle],'LineStyle',styles{i},'Color',colours{i}, 'DisplayName',['M_{NLE} ',names{i}]);
end
colormap(gca,R_colormap('Odds Ratio')); clim(Clims);
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 D(j).Ncb+1]); ylim([D(j).Mc,max([max(D(j).Ml2),max(D(j).Mlrg)]+0.3)]);
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvT plot.
ax2=subplot(2,4,7:8); hold on;
bp=bar([0,D(j).n]+1.5,D(j).Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(D(j).W)
    plot([1 D(j).Ncb+1],[i i]/length(D(j).W),'-w');
    bp(i).FaceColor=colours{i};
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1,D(j).Ncb+1]); ylim([0 1]);
set(gca, 'XScale', 'log');
linkaxes([ax1 ax2],'x');

% Plot earthquakes, stages, & well trajectory.
figure(1); clf;
% Map.
ax1=subplot(3,3,[1 2 4 5]); hold on;
plot(D(1).Wlon,   D(1).Wlat,'-k');
plot(D(1).Wlon(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    x=interp1(D(1).Wmd,D(1).Wlon,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(x,y,'-r','LineWidth',1);
    x=interp1(D(1).Wmd,D(1).Wlon,mean(S(i).Ploc,2),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,mean(S(i).Ploc,2),'linear');
    plot(x,y,'xr','LineWidth',1);
end
for i=1:length(D)
    scatter(D(i).Lon(D(i).Imb), D(i).Lat(D(i).Imb), D(i).Rs(D(i).Imb), D(i).OR(2:end),'filled','MarkerEdgeColor','k');
end
xlabel('Easting (m)'); ylabel('Northing (m)');
axis equal;
h = colorbar(); colormap(gca,R_colormap('Odds Ratio'));
clim(Clims);
ylabel(h, 'Log_{10} Odds Ratio');
% E-W depth profile.
ax2=subplot(3,3,[7 8]); hold on;
plot(D(1).Wlon,   D(1).Wtvd,'-k');
plot(D(1).Wlon(1),D(1).Wtvd(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    x=interp1(D(1).Wmd,D(1).Wlon,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(x,z,'-r','LineWidth',1);
    x=interp1(D(1).Wmd,D(1).Wlon,mean(S(i).Ploc,2),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,mean(S(i).Ploc,2),'linear');
    plot(x,z,'xr','LineWidth',1);
end
for i=1:length(D)
    scatter(D(i).Lon(D(i).Imb), -D(i).Dep(D(i).Imb), D(i).Rs(D(i).Imb), D(i).OR(2:end),'filled','MarkerEdgeColor','k');
end
xlabel('Easting (m)'); ylabel('Depth (m)');
set(gca, 'YDir','reverse');
colormap(gca,R_colormap('Odds Ratio'));
clim(Clims);
linkaxes([ax1 ax2],'x');
% N-S depth profile.
ax3=subplot(3,3,[3 6]); hold on;
plot(D(1).Wtvd,   D(1).Wlat,'-k');
plot(D(1).Wtvd(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    y=interp1(D(1).Wmd,D(1).Wlat,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(z,y,'-r','LineWidth',1);
    y=interp1(D(1).Wmd,D(1).Wlat,mean(S(i).Ploc,2),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,mean(S(i).Ploc,2),'linear');
    plot(z,y,'xr','LineWidth',1);
end
for i=1:length(D)
    scatter(-D(i).Dep(D(i).Imb), D(i).Lat(D(i).Imb), D(i).Rs(D(i).Imb), D(i).OR(2:end),'filled','MarkerEdgeColor','k');
end
xlabel('Depth (m)'); ylabel('Northing (m)');
colormap(gca,R_colormap('Odds Ratio'));
clim(Clims);
linkaxes([ax1 ax3],'y');






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end
