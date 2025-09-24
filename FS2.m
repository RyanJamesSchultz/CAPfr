% Script that visualizes the negative log-likelihood for Mmax fitting.
% Used to make Figure S2.
clear;
rng(7);

% Define some constants.
Nc=1e3;
Nr=1e0;
m1=0.0;
m2=2.4;
b=1.0;
d=0.01;

% Derive some variables.
m2=linspace(m2,m2,Nc);
a=log10(Nc);
Mmax=max(m2);

% Get a catalogue.
Mb=GR_MFD_Rand(m1,m2,  a,b, [1 Nc]);
Mu=GR_MFD_Rand(m1,Inf, a,b, [1 Nc]);
M2b=max(Mb):d:max(Mb)+2;
M2u=max(Mu):d:max(Mu)+2;

% Preallocate.
nLLb_gr=zeros(length(M2b),Nr);
nLLu_gr=nLLb_gr;
nLLb_dm=nLLb_gr;
nLLu_dm=nLLb_gr;

% Loop over all catalogue reshuffle trials.
for l=1:Nr
    l
    
    % Get a newly reshuffled catalogue.
    if(l==1)
        In=1:Nc;
    else
        In=randperm(Nc);
    end
    Mbi=Mb(In);
    Mui=Mu(In);
    
    % Collect the (bound) order statistics.
    [Mlrgb,Ib]=OrderStatistic(Mbi,Nc-0,'unique');
    Ib(isnan(Mlrgb))=[];
    Mlrgb(isnan(Mlrgb))=[];
    Mlrgb=[m1,Mlrgb];
    dMlrgb=diff(Mlrgb);
    Mlrgb(end)=[];

    % Collect the (unbound) order statistics.
    [Mlrgu,Iu]=OrderStatistic(Mui,Nc-0,'unique');
    Iu(isnan(Mlrgu))=[];
    Mlrgu(isnan(Mlrgu))=[];
    Mlrgu=[m1,Mlrgu];
    dMlrgu=diff(Mlrgu);
    Mlrgu(end)=[];

    % Loop over all of the possible Mmax values.
    for j=1:length(M2b)
        
        % Compute the GR-MFD negative log-likelihood.
        nLLb_gr(j,l)=-GR_MFD_LL(Mbi,m1,M2b(j),a,b)+GR_MFD_LL(Mbi,m1,Inf,a,b);
        nLLu_gr(j,l)=-GR_MFD_LL(Mui,m1,M2u(j),a,b)+GR_MFD_LL(Mui,m1,Inf,a,b);
        
        % Compute the dMlrg negative log-likelihood.
        nLLb_dm(j,l)=+GR_MFD_LL(dMlrgb,0,M2b(j)-Mlrgb,a,b)-GR_MFD_LL(dMlrgb,0,Inf,a,b);
        nLLu_dm(j,l)=+GR_MFD_LL(dMlrgu,0,M2u(j)-Mlrgu,a,b)-GR_MFD_LL(dMlrgu,0,Inf,a,b);
        
    end
end

% Get the composite negative log-likelihood.
nLLb=nLLb_gr+nLLb_dm;
nLLu=nLLu_gr+nLLu_dm;

% Get the MLE estimates of Mmax.
[nLLb_min,Ib]=min(nLLb,[],1);
[nLLu_min,Iu]=min(nLLu,[],1);




% Plot.

% Plots of the composite and GR-MFD nLL.
figure(52); clf;
% Plot the negative log-likelihood for a bound case.
subplot(211);
plot(M2b-M2b(1),nLLb,'-r'); hold on;
plot(M2b-M2b(1),nLLb_gr,'-b');
xlim([min(M2b) max(M2b)]-M2b(1));
Yl=ylim(); ylim([Yl(1) 1]);
plot(Mmax-M2b(1)*[1 1],ylim(),'-k');
plot(xlim(),[0 0],':k');
xlabel('Possible M_{MAX} Estimates (relative to M_{LRG})');
ylabel('Normalized Negative Log-Likelihood');
title('Bound Case');

% Plot the negative log-likelihood for an ubound case.
subplot(212);
plot(M2u-M2u(1),nLLu,'-r'); hold on;
plot(M2u-M2u(1),nLLu_gr,'-b');
xlim([min(M2u) max(M2u)]-M2u(1));
Yl=ylim(); ylim([Yl(1) -Yl(1)]);
plot(xlim(),[0 0],':k');
xlabel('Possible M_{MAX} Estimates (relative to M_{LRG})');
ylabel('Normalized Negative Log-Likelihood');
title('Unbound Case');




