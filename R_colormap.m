function [Z] = R_colormap(type_flag)
  % R_colormap()
  % 
  % Simple function to generate a color maps to display odds ratio information.
  %
  % Written by Ryan Schultz.
  %
  
  % Define the colormaps.
  if(strcmpi(type_flag,'Odds Ratio'))
      n=50;
      Z=[0       4     120;  ... % Dark blue.
         0       0     225;  ... % Blue.
         80      80    225;  ... % Light blue.
         255     255   255;  ... % White.
         255     255   255;  ... % White.
         255     255   255;  ... % Light magenta.
         255     176   225;  ... % Magenta.
         255     0     225;  ... % Magenta.
         150     0     150]; ... % Dark magenta.
     Z=flipud(Z);
  else
      n=50;
      Z=0;
  end
  
  % Make 16-bit RGB values between 0-1.
  Z=Z/255;
  
  % Interpolate to the desired number of samples.
  I=1:length(Z);
  Z=interp1(I,Z,min(I):1/n:max(I));
   
return;