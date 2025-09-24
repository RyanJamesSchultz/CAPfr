function [S]=EnsembleW(M,m1,S_m2,b,Nr,reshuffle_flag,smooth_flag)
  % Function to determine the sequence of ensemble weights for next largest 
  % events (NLEs), given the input catalogue and multiple model estimates 
  % of Mmax.
  % 
  % References:
  % Akaike (1998). Information theory and an extension of the maximum likelihood principle. In Selected papers of hirotugu akaike (pp. 199-213). New York, NY: Springer New York.
  % McQuarrie (1999). A small-sample correction for the Schwarz SIC model selection criterion. Statistics & probability letters, 44(1), 79-86, doi: 10.1016/S0167-7152(98)00294-6.
  % Schwarz (1978). Estimating the dimension of a model. The annals of statistics, 461-464.
  % Sugiura (1978). Further analysis of the data by akaike's information criterion and the finite corrections. Communications in Statistics-theory and Methods, 7(1), 13-26, doi: 10.1080/03610927808827599.
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Make M a row vector and get some important lengths.
  M=M(:)';
  Nc=length(M);
  Nm=length(S_m2);
  
  % Format the b-values appropriately.
  if(length(b)==length(M))
      b=b(:)';
  else
      b=b*ones(size(M));
  end
  
  % Predefine and preallocate the output structure and weight matrix.
  S=struct('Mmax',[],'W',[],'K',[]);
  for j=1:Nm
      if(length(S_m2(j).Mmax)==(Nc))
          S(j).Mmax=S_m2(j).Mmax(:)';
      else
          S(j).Mmax=S_m2(j).Mmax*ones([1 Nc]);
      end
      S(j).K=S_m2(j).K;
  end
  W=zeros([Nm Nc+1]);
  
  % Get the indicies of Mlrg.
  [~,I1]=OrderStatistic(M,Nc,'unique');
  
  % Loop over all reshuffling trials.
  for l=1:Nr+1
      
      % Get a newly reshuffled catalogue.
      if(l==1)
          Ml=M;
      else
          Il=randpermEW(Nc,I1,reshuffle_flag);
          Ml=M(Il);
      end
      
      % Compute the ensemble weights.
      Wl=getEW(Ml,m1,b,S);
      W=W+log(Wl)/(Nr+1);
  end
  W=exp(W);
  W=W./sum(W,1);
  
  % Optionally, smooth out the sequence of weights temporally.
  wp=W(:,1);
  for i=2:Nc
      w=W(:,i);
      if(strcmpi(smooth_flag,'mean'))
          w=w+wp;
          w=w/sum(w);
      elseif(strcmpi(smooth_flag,'geomean'))
          w=sqrt(w.*wp);
          w=w/sum(w);
      elseif(strcmpi(smooth_flag,'both'))
          w1=w+wp;        w1=w1/sum(w1);
          w2=sqrt(w.*wp); w2=w2/sum(w2);
          w=w1+w2;
          w=w/sum(w);
      end
      wp=w;
      W(:,i)=w;
  end
  
  % Stuff the results into the output data structure.
  for j=1:Nm
      S(j).W=W(j,:);
  end
  
end




%%%% SUBROUNTINES.

% The ensemble weights for dMlrg(i).
function [W]=getEW(M,m1,b,S)
  
  % Get the Mlrg & dMlrg sequences (and indicies to when they occur).
  [Mlrg,I]=OrderStatistic(M,length(M)-0,'unique');
  I(isnan(Mlrg))=[];
  Mlrg(isnan(Mlrg))=[];
  Mlrg=[m1,Mlrg];
  dMlrg=diff(Mlrg);
  Mlrg(end)=[];
  
  % Define relevant lengths.
  Nc=length(M);
  Nm=length(S);
  
  % Predefine the log-likelihood vector.
  LL=zeros([Nm Nc]);
  
  % Loop over all of the proposed Mmax models.
  for k=1:Nm
      
      % Get the proposed Mmax.
      m2=S(k).Mmax;
      
      % Construct the composite log-likelihood.
      LL(k,:)=log(GR_MFD(M,m1,m2,0,b,'norm')); % GR-MFD term.
      LL(k,I)=LL(k,I)-log(GR_MFD(dMlrg,0,m2(I)-Mlrg,0,b(I),'norm')); % dMlrg term (approximate).
      
  end
  LL=cumsum(LL,2);
  
  % Make parameter & sample size into apprpriately sized matrices. 
  %n=1:Nc;
  %K=3+1; % b,m1,m2.
  n=repmat(1:Nc,  [Nm 1]);
  K=repmat([S.K]',[1 Nc]);
  
  % Get the AIC and BIC scores.
  %AIC=2*K-2*LL; % [Akaike, 1998].
  %BIC=log(n).*K-2*LL; % [Schwarz, 1978].
  AICc=2*K+(2*K.*(K+1)./(n-K-1))-2*LL; % [Sugiura, 1978].
  BICc=n.*log(n).*K./(n-K-1)-2*LL; % [McQuarrie, 1999].
  
  % Do a smooth/linear interpolation for when the sample size is too small.
  I=(n<=K+2);
  AICc(I)=((2*K(I)+(2*K(I).*(K(I)+1)))./(K(I)+2)).*n(I)-2*LL(I);
  BICc(I)=(K(I).*log(K(I)+2)).*n(I)-2*LL(I);
  
  % Subtract off the best model's score.
  dAIC=AICc-min(AICc,[],1);
  dBIC=BICc-min(BICc,[],1);
  
  % Get the ensemble weights [Wagenmakers & Farrell, 2004].
  Wa=exp(-dAIC/2);
  Wb=exp(-dBIC/2);
  Wa=Wa./sum(Wa,1);
  Wb=Wb./sum(Wb,1);
  
  % Average between the two weights.
  W=exp(log(Wa)+log(Wb));
  %W=Wa;
  %W=Wb;
  W=W./sum(W,1);
  W=[ones([Nm 1])/Nm, W];
  
end


% Makes index permutation, with Mlrg constraints.
function [Ir]=randpermEW(n,I,reshuffle_flag)
  
  % Define the index blocks and predefine the output vector.
  I=I(:)';
  dI=diff([I,n+1]);
  Ir=[];
  
  % Simply return a regular permutation, if flagged to.
  if(strcmpi(reshuffle_flag,'regular'))
      Ir=randperm(n);
      return;
  end
  
  % Loop over all blocks, and then permute within each of them.
  for i=1:length(dI)
      if(strcmpi(reshuffle_flag,'blocked'))
          Ii=randperm(dI(i))+I(i)-1; % The whole block permutes.
          Ir=[Ir,Ii];
      elseif(strcmpi(reshuffle_flag,'fixed'))
          Ii=randperm(dI(i)-1)+I(i); % Mlrg stays fixed, the rest of the block permutes.
          Ir=[Ir,I(i),Ii];
      end
  end
  
end

