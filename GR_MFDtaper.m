function [PDF,CDF,SvF]=GR_MFDtaper(M,m1,m2,a,b,norm_flag)
  % Function that computes the tapered Gutenberg-Richter magnitude-frequency 
  % distributions (GR-MFD).
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Convert magnitudes to moment.
  dm2=0.45;
  Mo=10.^(1.5*M+9.1); % Mw to Mo (Nm).
  M1=10.^(1.5*m1+9.1);
  M2=10.^(1.5*(m2)+9.1);
  %M2b=10.^(1.5*(m2)+9.1);
  
  % Define some useful values.
  beta=2*b/3;
  
  % The PDF & CDF for a doubly bounded GR-MFD.
  SvF=exp((M1-Mo)./M2).*(M1./Mo).^beta;   SvF(Mo<M1)=NaN; %SvF(Mo>M2b)=NaN;
  PDF=SvF.*((beta./Mo)+(1./M2));          PDF(Mo<M1)=NaN; %PDF(Mo>M2b)=NaN;
  CDF=(1-SvF);                            CDF(Mo<M1)=NaN; %CDF(Mo>M2b)=NaN;
  
  % Change the normalization, if flagged to.
  if(strcmpi(norm_flag,'count'))
      PDF=PDF.*10.^a;
      CDF=CDF.*10.^a;
      SvF=SvF.*10.^a;
  end
  
end