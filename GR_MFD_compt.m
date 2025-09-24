function [PDF,CDF,SvF]=GR_MFD_compt(M,m1,m2,a,b)
  % Function that computes the effective Gutenberg-Richter 
  % magnitude-frequency distribution (GR-MFD), if both/either the b-value 
  % or Mmax are temporally varying.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the upper/lower magnitude bounds.
  PDF=zeros(size(M)); CDF=PDF; SvF=PDF;
  Nm=length(m2);
  
  % Loop over all of the event data and .
  for i=1:length(m2)
      [pdf,cdf,svf]=GR_MFDtaper(M,m1,m2(i),1,b(i),'norm');
      pdf(isnan(pdf))=0; cdf(isnan(cdf))=1; svf(isnan(svf))=0;
      PDF=PDF+pdf/Nm; CDF=CDF+cdf/Nm; SvF=SvF+svf/Nm;
  end
  
end