% Locations and timing of earthquakes and stages at Helsinki St1.
% Used to make Figure S11.
clear;

% Get all of the case data.
D= PreProcData({'St1'},1:5); % Cluster 1.

% Stage locations.
%D(1).Slon=interp1(D(1).Wmd,D(1).Wlon,D(1).MDs,'linear'); D(1).Slat=interp1(D(1).Wmd,D(1).Wlat,D(1).MDs,'linear'); D(1).Stvd=interp1(D(1).Wmd,D(1).Wtvd,D(1).MDs,'linear');

% Define some colours and plotting details that I'd like to use.
colours={'#d8973c'};
ColorInj='#87b6e1';
D(1).Rs=getMscale(D(1).M)*1.0;

% Make well surface location the origin.
%Xc=D(1).Wlon(1);
%Yc=D(1).Wlat(1);
Zc=D(1).Wtvd(1);
Xc=2.54904e7;
Yc=6.67505e6;




% Plot.
figure(511); clf;

% MvT plot (full).
subplot(311); hold on;
area(D(1).t,D(1).v*5-1,-1,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2018,06,03,00,00,00) datetime(2018,08,01,12,00,00)]); ylim([-1 2]);

% Map.
subplot(3,1,[2 3]); hold on;
scatter(D(1).Lon-Xc,D(1).Lat-Yc,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
%plot(D(1).Slon-Xc,D(1).Slat-Yc,'dk','MarkerFaceColor',colours{1});
% Well PNR-1z.
plot(D(1).Wlon-Xc,   D(1).Wlat-Yc,'-k');
plot(D(1).Wlon(1)-Xc,D(1).Wlat(1)-Yc,'ok','MarkerFaceColor','k');
xlabel('Easting (m)'); ylabel('Northing (m)');
xlim([-250 1050]); ylim([-200 1100]);





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end