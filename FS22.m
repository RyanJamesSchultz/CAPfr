% Locations and timing of earthquakes and stages at PNR-2.
% Used to make Figure S22.
clear;

% Get all of the case data.
D1=PreProcData({'PNR2'},1:3); % Cluster W.
D2=PreProcData({'PNR2'},4);   % Stage 4.
D3=PreProcData({'PNR2'},5:7); % Cluster E.
D(1)=D1; D(2)=D2; D(3)=D3;

% Stage locations.
D(1).Slon=interp1(D(1).Wmd,D(1).Wlon,D(1).MDs,'linear'); D(1).Slat=interp1(D(1).Wmd,D(1).Wlat,D(1).MDs,'linear'); D(1).Stvd=interp1(D(1).Wmd,D(1).Wtvd,D(1).MDs,'linear');
D(2).Slon=interp1(D(2).Wmd,D(2).Wlon,D(2).MDs,'linear'); D(2).Slat=interp1(D(2).Wmd,D(2).Wlat,D(2).MDs,'linear'); D(2).Stvd=interp1(D(2).Wmd,D(2).Wtvd,D(2).MDs,'linear');
D(3).Slon=interp1(D(3).Wmd,D(3).Wlon,D(3).MDs,'linear'); D(3).Slat=interp1(D(3).Wmd,D(3).Wlat,D(3).MDs,'linear'); D(3).Stvd=interp1(D(3).Wmd,D(3).Wtvd,D(3).MDs,'linear');

% Define some colours and plotting details that I'd like to use.
colours={'#d8973c', '#6c6d6c', '#a1cda8'};
ColorInj='#87b6e1';
D(1).Rs=getMscale(D(1).M)*1.0;
D(2).Rs=getMscale(D(2).M)*1.0;
D(3).Rs=getMscale(D(3).M)*1.0;



% Plot.
figure(522); clf;

% Use well surface location as the figure origin.
Xc=D(1).Wlon(1); Yc=D(1).Wlat(1); Zc=D(1).Wtvd(1); 

% MvT plot (2022).
subplot(311); hold on;
% CW.
area(D(1).t,D(1).v/1-2,-2,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
% S4.
area(D(2).t,D(2).v/1-2,-2,'FaceColor',ColorInj);
scatter(D(2).T,D(2).M,D(2).Rs,'MarkerFaceColor',colours{2},'MarkerEdgeColor','none');
plot(D(2).T,cummax(D(2).M),'-','Color',colours{2});
% CE.
area(D(3).t,D(3).v/1-2,-2,'FaceColor',ColorInj);
scatter(D(3).T,D(3).M,D(3).Rs,'MarkerFaceColor',colours{3},'MarkerEdgeColor','none');
plot(D(3).T,cummax(D(3).M),'-','Color',colours{3});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2019,08,15,00,00,00) datetime(2019,08,28,00,00,00)]); ylim([-2 3]);

% Map.
subplot(3,1,[2 3]); hold on;
% CW.
scatter(D(3).Lon-Xc,D(3).Lat-Yc,D(3).Rs,'MarkerFaceColor',colours{3},'MarkerEdgeColor','none');
plot(D(3).Slon-Xc,D(3).Slat-Yc,'dk','MarkerFaceColor',colours{3});
% S4.
scatter(D(2).Lon-Xc,D(2).Lat-Yc,D(2).Rs,'MarkerFaceColor',colours{2},'MarkerEdgeColor','none');
plot(D(2).Slon-Xc,D(2).Slat-Yc,'dk','MarkerFaceColor',colours{2});
% CE.
scatter(D(1).Lon-Xc,D(1).Lat-Yc,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).Slon-Xc,D(1).Slat-Yc,'dk','MarkerFaceColor',colours{1});
% Well PNR-2.
plot(D(1).Wlon-Xc,   D(1).Wlat-Yc,'-k');
plot(D(1).Wlon(1)-Xc,D(1).Wlat(1)-Yc,'ok','MarkerFaceColor','k');
xlabel('Easting (m)'); ylabel('Northing (m)');





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end