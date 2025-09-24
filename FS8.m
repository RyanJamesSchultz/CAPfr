% Locations and timing of earthquakes and stages at PNR-1z.
% Used to make Figure S8.
clear;

% Get all of the case data.
D= PreProcData({'PNR1z'},1:23); % Cluster 1.

% Stage locations.
D(1).Slon=interp1(D(1).Wmd,D(1).Wlon,D(1).MDs,'linear'); D(1).Slat=interp1(D(1).Wmd,D(1).Wlat,D(1).MDs,'linear'); D(1).Stvd=interp1(D(1).Wmd,D(1).Wtvd,D(1).MDs,'linear');

% Define some colours and plotting details that I'd like to use.
colours={'#d8973c'};
ColorInj='#87b6e1';
D(1).Rs=getMscale(D(1).M)*1.0;

% Make well surface location the origin.
Xc=D(1).Wlon(1);
Yc=D(1).Wlat(1);
Zc=D(1).Wtvd(1);





% Plot.
figure(58); clf;

% MvT plot (full).
subplot(511); hold on;
area(D(1).t,D(1).v/3-2,-2,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2018,10,16,00,00,00) datetime(2018,12,18,12,00,00)]); ylim([-2 2]);

% MvT plot (1st half).
subplot(512); hold on;
area(D(1).t,D(1).v/3-2,-2,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2018,10,16,00,00,00) datetime(2018,11,04,00,00,00)]); ylim([-2 2]);

% MvT plot (2nd half).
subplot(513); hold on;
area(D(1).t,D(1).v/3-2,-2,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2018,12,08,00,00,00) datetime(2018,12,18,12,00,00)]); ylim([-2 2]);

% Map.
subplot(5,1,[4 5]); hold on;
scatter(D(1).Lon-Xc,D(1).Lat-Yc,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).Slon-Xc,D(1).Slat-Yc,'dk','MarkerFaceColor',colours{1});
% Well PNR-1z.
plot(D(1).Wlon-Xc,   D(1).Wlat-Yc,'-k');
plot(D(1).Wlon(1)-Xc,D(1).Wlat(1)-Yc,'ok','MarkerFaceColor','k');
xlabel('Easting (m)'); ylabel('Northing (m)');
xlim([-1900 100]); ylim([-300 200]);





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end