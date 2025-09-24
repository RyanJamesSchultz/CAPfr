% Simple script to visualize the perf/stage spacing at FORGE well 16A.
% Used to make Figure S14.
clear;

% Define an input structure.
S=struct('Stage',[],'Ploc',[]);

% Define the perf/stage locations.
S( 1).Stage=1;  S( 1).Ploc(1,:)=[10787 10955]*0.3048;
S( 2).Stage=2;  S( 2).Ploc(1,:)=[10560 10580]*0.3048;
S( 3).Stage=3;  S( 3).Ploc(1,:)=[10120 10140]*0.3048;
S( 4).Stage=4;  S( 4).Ploc(1,:)=[10070 10076]*0.3048;
S( 5).Stage=5;  S( 5).Ploc(1,:)=[10020 10026]*0.3048;
S( 6).Stage=6;  S( 6).Ploc(1,:)=[ 9970  9976]*0.3048; S( 6).Ploc(2,:)=[9959 9962]*0.3048;
S( 7).Stage=7;  S( 7).Ploc(1,:)=[ 9898  9901]*0.3048; S( 7).Ploc(2,:)=[9850 9853]*0.3048; S( 7).Ploc(3,:)=[9798 9801]*0.3048;
S( 8).Stage=8;  S( 8).Ploc(1,:)=[ 9720  9723]*0.3048; S( 8).Ploc(2,:)=[9695 9698]*0.3048; S( 8).Ploc(3,:)=[9670 9673]*0.3048; S( 8).Ploc(4,:)=[9645 9648]*0.3048; S( 8).Ploc(5,:)=[9620 9623]*0.3048; S( 8).Ploc(6,:)=[9595 9598]*0.3048; S( 8).Ploc(7,:)=[9570 9573]*0.3048; S( 8).Ploc(8,:)=[9545 9548]*0.3048;
S( 9).Stage=9;  S( 9).Ploc(1,:)=[ 9490  9493]*0.3048; S( 9).Ploc(2,:)=[9470 9473]*0.3048; S( 9).Ploc(3,:)=[9445 9448]*0.3048; S( 9).Ploc(4,:)=[9420 9423]*0.3048; S( 9).Ploc(5,:)=[9395 9398]*0.3048; S( 9).Ploc(6,:)=[9370 9373]*0.3048; S( 9).Ploc(7,:)=[9345 9348]*0.3048; S( 9).Ploc(8,:)=[9320 9323]*0.3048;
S(10).Stage=10; S(10).Ploc(1,:)=[ 9270  9276]*0.3048;

% Load in the well data.
D=PreProcData({'FORGE24'},1);


% Plot well trajectory.
figure(514); clf;
% Map.
ax1=subplot(3,3,[1 2 4 5]);
plot(D(1).Wlon,   D(1).Wlat,'-k'); hold on;
plot(D(1).Wlon(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    x=interp1(D(1).Wmd,D(1).Wlon,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(x,y,'-r','LineWidth',1);
    x=interp1(D(1).Wmd,D(1).Wlon,mean(S(i).Ploc,2),'linear');
    y=interp1(D(1).Wmd,D(1).Wlat,mean(S(i).Ploc,2),'linear');
    plot(x,y,'xr','LineWidth',1);
end
xlabel('Easting (m)'); ylabel('Northing (m)');
% E-W depth profile.
ax2=subplot(3,3,[7 8]);
plot(D(1).Wlon,   D(1).Wtvd,'-k'); hold on;
plot(D(1).Wlon(1),D(1).Wtvd(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    x=interp1(D(1).Wmd,D(1).Wlon,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(x,z,'-r','LineWidth',1);
    x=interp1(D(1).Wmd,D(1).Wlon,mean(S(i).Ploc,2),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,mean(S(i).Ploc,2),'linear');
    plot(x,z,'xr','LineWidth',1);
end
xlabel('Easting (m)'); ylabel('Depth (m)');
set(gca, 'YDir','reverse');
linkaxes([ax1 ax2],'x');
% N-S depth profile.
ax3=subplot(3,3,[3 6]);
plot(D(1).Wtvd,   D(1).Wlat,'-k'); hold on;
plot(D(1).Wtvd(1),D(1).Wlat(1),'ok','MarkerFaceColor','k');
for i=1:length(S)
    y=interp1(D(1).Wmd,D(1).Wlat,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,min(S(i).Ploc(:)):max(S(i).Ploc(:)),'linear');
    plot(z,y,'-r','LineWidth',1);
    y=interp1(D(1).Wmd,D(1).Wlat,mean(S(i).Ploc,2),'linear');
    z=interp1(D(1).Wmd,D(1).Wtvd,mean(S(i).Ploc,2),'linear');
    plot(z,y,'xr','LineWidth',1);
end
xlabel('Depth (m)'); ylabel('Northing (m)');
linkaxes([ax1 ax3],'y');




% Determine some values of interest.
mp=[];
for i=1:length(S)
    ms(i)=mean(S(i).Ploc(:));
    mp=[mp;mean(S(i).Ploc,2)]
end
