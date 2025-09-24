% Show practical implications of CAP-test results.
% Used to make Figures 7 & S42.
clear;

% Predefine some constants.
thresh=3;
Mc=0.0;
dMg=0.10;
dMr=0.30;

% Make list of cases/clusters to consider.
CaseList={};
CaseList=[CaseList, {'FORGE-c1.mat','FORGE-c2a.mat','FORGE-c2.mat','FORGE-c3.mat'}];
CaseList=[CaseList, {'PNR1z-all.mat','PNR2-cE.mat','PNR2-cW.mat','PNR2-s4.mat'}];
CaseList=[CaseList, {'SSFS93-all.mat','SSFS00-all.mat','SSFS03-all.mat','SSFS04-all.mat','SSFS05-all.mat'}];
CaseList=[CaseList, {'St1-all.mat'}];

% Preallocate.
dMb=[];
dMu=[];

% Loop over all of the cases.
figure(542); clf;
for i=1:length(CaseList)
    
    % Load data and then get the catalogue NLE jumps.
    load(CaseList{i},'D');
    [dMlrg1,dMlrg2,j,mOR]=get_dMconf(D,thresh);
    
    % Stuff the results into the correct vector.
    if(isempty(j))
        dMu=[dMu;dMlrg1];
    else
        dMu=[dMu;dMlrg1];
        dMb=[dMb;dMlrg2];
    end
    
    % Add to the plot.
    Name=regexp(CaseList{i}, '[-.]', 'split');
    l=loglog(0:(length(mOR)-1), mOR, 'DisplayName',[Name{1},' ',Name{2}]); hold on;
    loglog(j-1, mOR(j),'o','MarkerEdgeColor',l.Color, 'MarkerFaceColor',l.Color, 'HandleVisibility','off');
    
end

% Clean up plot.
loglog(xlim(), thresh*[1 1], '--k', 'HandleVisibility','off');
loglog(xlim(), 1e0*[1 1], '-k', 'HandleVisibility','off');
loglog(xlim(), 3e0*[1 1], ':k', 'HandleVisibility','off');
loglog(xlim(), 1e1*[1 1], ':k', 'HandleVisibility','off');
loglog(xlim(), 1e2*[1 1], ':k', 'HandleVisibility','off');
ylim([1e-3 1e+3]);
xlabel('Chronoloigcal Event Count');
ylabel('Relative Odds Ratio');
legend('Location','northwest');

% Make a synthetic catalogue.
[dMs]=GR_MFD_Rand(Mc,Inf,1,1.0,size(dMb));

% Do some KS-tests.
[~,pv_bu]=kstest2(dMb,dMu);
[~,pv_bs]=kstest2(dMb,dMs);
[~,pv_us]=kstest2(dMu,dMs);

% Fit the GR-MFDs.
[bb,bb_err,ab,R2b,~,Mgrb,Ngrb,ngrb]=Bval(dMb,Mc,dMg);
[bu,bu_err,au,R2u,~,Mgru,Ngru,ngru]=Bval(dMu,Mc,dMg);
[bs,bs_err,as,R2s,~,Mgrs,Ngrs,ngrs]=Bval(dMs,Mc,dMg);

% Get the GR-MFD plotting values.
pb=[-mean(bb),ab];
Mgr_fitb=[Mc, max(dMb)];
Ngr_fitb=10.^polyval(pb,Mgr_fitb);
pu=[-mean(bu),au];
Mgr_fitu=[Mc, max(dMu)];
Ngr_fitu=10.^polyval(pu,Mgr_fitu);
ps=[-mean(bs),as];
Mgr_fits=[Mc, max(dMs)];
Ngr_fits=10.^polyval(ps,Mgr_fits);

% Report some values.
1-pv_bu
[bb bb_err]
R2b
[bu bu_err]
R2u
median(dMb)




% Plot.

% Plot the dMLRG histograms, in the style of GR-MFDs.
GREY=[0.85,0.85,0.85];
figure(7); clf;
% Bound and unbound data.
%subplot(211);
semilogy(Mgrb, Ngrb, 'o', 'Color', 'b'); hold on;
bar(Mgrb,ngrb, 'FaceColor', 'b');
semilogy(Mgr_fitb, Ngr_fitb, '-', 'Color', 'b');
semilogy(Mgru, Ngru, 'o', 'Color', 'r'); hold on;
bar(Mgru,ngru, 'FaceColor', 'r');
semilogy(Mgr_fitu, Ngr_fitu, '-', 'Color', 'r');
xlim([min(Mgrb)-dMg/2 max([Mgrb,Mgru])+dMg/2]); ylim([0.8 1.3*max([Ngrb,Ngru])]);
plot(dMr*[1 1],ylim,'--k');
xlabel('\DeltaM_{LRG}'); ylabel('Count');
% Synthetic data.
%subplot(212);
%semilogy(Mgrs, Ngrs, 'o', 'Color', 'k'); hold on;
%bar(Mgrs,ngrs, 'FaceColor', GREY);
%semilogy(Mgr_fits, Ngr_fits, '-', 'Color', 'black');
%xlim([min(Mgrs)-dMg/2 max([Mgrb,Mgru,Mgrs])+dMg/2]); ylim([0.7 1.3*max(Ngrs)]);
%plot(0*[1 1],ylim,'--k');
%xlabel('\DeltaM_{LRG}'); ylabel('Count');

% Get the confusion matrix.
Cm=ConfusionMatrix(dMb,dMu,dMr);

% Plot.
figure(8); clf;
cc=confusionchart(Cm,{'Bound','Unbound'});
%cc.RowSummary = 'row-normalized';
cc.ColumnSummary = 'column-normalized';
xlabel('TLP Result'); ylabel('EW-test Result');






%%%% SUBROUNTINES.

% Get the jumps in magnitude, after the EW-test hits an input threshold.
function [dMlrg1,dMlrg2,i,tOR]=get_dMconf(D,threshold)
  
  % Get the catalogue observables.
  M=D.M(D.M>=D.Mc);
  Mlrg=OrderStatistic(M,length(M)-0,'none'); Mlrg=[D.Mc;Mlrg];
  
  % Find the event index that crosses the threshold.
  tOR=max([D.Wb(1,:)./D.Wb(3,:); D.Wb(2,:)./D.Wb(3,:)]);
  i=find(tOR>=threshold,1,'first');
  
  % Split the catalogue.
  if(isempty(i))
      dMlrg1=diff(unique(Mlrg));
      dMlrg2=[];
  else
      dMlrg1=diff(unique(Mlrg(1:i)));
      dMlrg2=diff(unique(Mlrg(i:end)));
  end
  
end

% Make a confusion matrix for cases.
function [Cm]=ConfusionMatrix(dMb,dMu,threshold)
  
  % Populate the matrix.
  Cm=zeros(2);
  Cm(1,1)=sum(dMb<threshold);
  Cm(2,1)=sum(dMb>=threshold);
  Cm(1,2)=sum(dMu<threshold);
  Cm(2,2)=sum(dMu>=threshold);
  
end