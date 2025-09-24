% Locations and timing of earthquakes and stages at FORGE.
% Used to make Figure S15.
clear;

% Get all of the case data.
D1 =PreProcData({'FORGE22'},1:2); % Cluster 1.
D2a=PreProcData({'FORGE22'},3);   % Cluster 2a.
D2 =PreProcData({'FORGE24'},1:4); % Cluster 2.
D3 =PreProcData({'FORGE24'},5:8); % Cluster 3.
%D4= PreProcData({'FORGE24'},5:8); % Cluster 4.
%D5= PreProcData({'FORGE24'},4);   % Stage 6.
D(1)=D1; D(2)=D2a; D(3)=D2; D(4)=D3; %D(5)=D4; D(6)=D5;

% Stage locations.
D(1).Slon=interp1(D(1).Wmd,D(1).Wlon,D(1).MDs,'linear'); D(1).Slat=interp1(D(1).Wmd,D(1).Wlat,D(1).MDs,'linear'); D(1).Stvd=interp1(D(1).Wmd,D(1).Wtvd,D(1).MDs,'linear');
D(2).Slon=interp1(D(2).Wmd,D(2).Wlon,D(2).MDs,'linear'); D(2).Slat=interp1(D(2).Wmd,D(2).Wlat,D(2).MDs,'linear'); D(2).Stvd=interp1(D(2).Wmd,D(2).Wtvd,D(2).MDs,'linear');
D(3).Slon=interp1(D(3).Wmd,D(3).Wlon,D(3).MDs,'linear'); D(3).Slat=interp1(D(3).Wmd,D(3).Wlat,D(3).MDs,'linear'); D(3).Stvd=interp1(D(3).Wmd,D(3).Wtvd,D(3).MDs,'linear');
D(4).Slon=interp1(D(4).Wmd,D(4).Wlon,D(4).MDs,'linear'); D(4).Slat=interp1(D(4).Wmd,D(4).Wlat,D(4).MDs,'linear'); D(4).Stvd=interp1(D(4).Wmd,D(4).Wtvd,D(4).MDs,'linear');
%D(5).Slon=interp1(D(5).Wmd,D(5).Wlon,D(5).MDs,'linear'); D(5).Slat=interp1(D(5).Wmd,D(5).Wlat,D(5).MDs,'linear'); D(5).Stvd=interp1(D(5).Wmd,D(5).Wtvd,D(5).MDs,'linear');
%D(6).Slon=interp1(D(6).Wmd,D(6).Wlon,D(6).MDs,'linear'); D(6).Slat=interp1(D(6).Wmd,D(6).Wlat,D(6).MDs,'linear'); D(6).Stvd=interp1(D(6).Wmd,D(6).Wtvd,D(6).MDs,'linear');

% Define some colours and plotting details that I'd like to use.
colours={'#d8973c', '#bed3da', '#628395', '#a1cda8', '#9000ff', '#6c6d6c'};
ColorInj='#87b6e1';
D(1).Rs=getMscale(D(1).M)*1.0;
D(2).Rs=getMscale(D(2).M)*1.0;
D(3).Rs=getMscale(D(3).M)*1.0;
D(4).Rs=getMscale(D(4).M)*1.0;
%D(5).Rs=getMscale(D(5).M)*1.0;
%D(6).Rs=getMscale(D(6).M)*1.0;



% Plot.
figure(515); clf;

% MvT plot (2022).
subplot(411); hold on;
% C1.
area(D(1).t,D(1).v/3-2,-2,'FaceColor',ColorInj);
scatter(D(1).T,D(1).M,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).T,cummax(D(1).M),'-','Color',colours{1});
% C2a.
area(D(2).t,D(2).v/3-2,-2,'FaceColor',ColorInj);
scatter(D(2).T,D(2).M,D(2).Rs,'MarkerFaceColor',colours{2},'MarkerEdgeColor','none');
plot(D(2).T,cummax(D(2).M),'-','Color',colours{2});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2022,04,17,00,00,00) datetime(2022,04,23,00,00,00)]); ylim([-2 1]);

% MvT plot (2024).
subplot(412); hold on;
% C2.
area(D(3).t,D(3).v/6-1/2,-1/2,'FaceColor',ColorInj);
scatter(D(3).T,D(3).M,D(3).Rs,'MarkerFaceColor',colours{3},'MarkerEdgeColor','none');
plot(D(3).T,cummax(D(3).M),'-','Color',colours{3});
% C3.
area(D(4).t,D(4).v/6-1/2,-1/2,'FaceColor',ColorInj);
scatter(D(4).T,D(4).M,D(4).Rs,'MarkerFaceColor',colours{4},'MarkerEdgeColor','none');
plot(D(4).T,cummax(D(4).M),'-','Color',colours{4});
% S6.
%area(D(6).t,D(6).v/5-1,-1,'FaceColor',ColorInj);
%scatter(D(6).T,D(6).M,D(6).Rs,'MarkerFaceColor',colours{6},'MarkerEdgeColor','none');
%plot(D(6).T,cummax(D(6).M),'-','Color',colours{6});
% C4.
%area(D(5).t,D(5).v/5-1,-1,'FaceColor',ColorInj);
%scatter(D(5).T,D(5).M,D(5).Rs,'MarkerFaceColor',colours{5},'MarkerEdgeColor','none');
%plot(D(5).T,cummax(D(5).M),'-','Color',colours{5});
xlabel('Time'); ylabel('Magnitude');
xlim([datetime(2024,04,03,12,00,00) datetime(2024,04,07,16,00,00)]); ylim([-1/2 2]);

% Map.
subplot(4,1,[3 4]); hold on;
% C4.
%scatter(D(5).Lon,D(5).Lat,D(5).Rs,'MarkerFaceColor',colours{5},'MarkerEdgeColor','none');
%plot(D(5).Slon,D(5).Slat,'dk','MarkerFaceColor',colours{5});
% S6.
%scatter(D(6).Lon,D(6).Lat,D(6).Rs,'MarkerFaceColor',colours{6},'MarkerEdgeColor','none');
%plot(D(6).Slon,D(6).Slat,'dk','MarkerFaceColor',colours{6});
% C3.
scatter(D(4).Lon,D(4).Lat,D(4).Rs,'MarkerFaceColor',colours{4},'MarkerEdgeColor','none');
plot(D(4).Slon,D(4).Slat,'dk','MarkerFaceColor',colours{4});
% C3a.
scatter(D(3).Lon,D(3).Lat,D(3).Rs,'MarkerFaceColor',colours{3},'MarkerEdgeColor','none');
plot(D(3).Slon,D(3).Slat,'dk','MarkerFaceColor',colours{3});
% C2.
scatter(D(2).Lon,D(2).Lat,D(2).Rs,'MarkerFaceColor',colours{2},'MarkerEdgeColor','none');
plot(D(2).Slon,D(2).Slat,'dk','MarkerFaceColor',colours{2});
% C1.
scatter(D(1).Lon,D(1).Lat,D(1).Rs,'MarkerFaceColor',colours{1},'MarkerEdgeColor','none');
plot(D(1).Slon,D(1).Slat,'dk','MarkerFaceColor',colours{1});
% Well 16A & 16B
plot(D(1).Wlon,   D(1).Wlat,'-k');
plot(D(1).Wlon(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
xlabel('Easting (m)'); ylabel('Northing (m)');





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end