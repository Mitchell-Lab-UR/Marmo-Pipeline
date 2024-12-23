function [tcor,acor,pcor,ncor,BI,PK,LENO] = comp_autocor_fit(sptimes,binsize,maxlag)
% function [tcor,acor,pcor,ncor,BI,PK,LENO] = comp_autocor_fit(sptimes,binsize,maxlag)
%
%** inputs:  sptimes, Nx1 list of spike times in secs
%                if empty, it will build an artificial spike train
%                    with a refractory and burst process
%            binsize, 0.2 ms default, size of bins for spikes
%            maxlag, maximum lag to show autocor out to, 30ms default
%*** outputs:  returns autocorrelation, tcor is timelags (0 to maxlag), 
%***                                    acor is autocorrelation 
%***                                    pcor is poisson expectation of auto
%***                                    ncor - smoothed acor
%***                                    BI - computed burst index
%***                                    PK - peak time of smoothed autocor
%***                                    LENO - number of spikes used
%***
%***  This version will fit a log-gauss basis set to the autocorr,
%***  which should help for units of low firing rates,
%***  and then use that to get BI and PK values

  MAXSPIKES = 50000;  % do not process more than this (takes too long)
  plotprogress = 0;   % show results as you go
                                            
  if isempty(binsize)
      BinSize = 0.0004;  % 0.4 ms (in secs)
  else
      BinSize = binsize;
  end 
  if isempty(maxlag)
      MaxLag = 0.040;    % 40 ms lag (in secs)
  else
      MaxLag = maxlag;
  end
  
  acor = 1;
  tcor = 1;
  ncor = 1;
  pcor = 1;
 
  %******* use xcorr to compute autocorrelation on binned spikes
  disp('Computing correlation over lags');
  NLag = 1+floor(MaxLag/BinSize);
  Asum = zeros(1,NLag);
  LENO = min(length(sptimes),MAXSPIKES);
  sptimes = sptimes(1:LENO);  % truncate to reduce search speed in computing
  k10 = floor(LENO/10);  % report progress in 10ths
  for k = 1:LENO
      spt = sptimes(k);
      %** any time zero multiplies, you get zero
      %** thus you only need those moments that bin 0 has a one
      zz = find( (sptimes > spt) & (sptimes < (spt+MaxLag+BinSize)) );
      if ~isempty(zz)
          tlags = 1 + floor( (sptimes(zz)-spt)/BinSize );
          Asum(tlags) = Asum(tlags)+1;
      end
      %***  report progress on command line
      if (mod(k,k10) == 0)
          disp(sprintf('Progress %d percent',ceil(k*100/LENO)));
      end
      %*******
  end
  acor = Asum / LENO; 
  tcor = (1000 * BinSize) * (0:(NLag-1));
  %****** normalize by expectation for a Poisson process (flat autocor)
  disp('Computing Poisson expectation');
  GoTime = min(sptimes);
  FiTime = max(sptimes);
  rate = length(sptimes)/(FiTime-GoTime);
  arate = (rate*BinSize);
  pcor = arate * ones(size(acor));
  % ncor = acor ./ pcor;
  %*** Gaussian smooth autocor function
  % ncor = gauss_smooth(ncor,floor(0.0012/BinSize));
  yy = basis(8,MaxLag,BinSize);
  b = regress(acor',yy');  %using least-squares ... better if log like?
  ncor = (yy' * b)';
  ncor = ncor ./ pcor;
  z = find( ncor == max(ncor));
  PK = tcor(z(1));  % peak time of autocor
  BI = max(ncor);
  
  %********
  % z1 = find( (tcor >= 1) & (tcor < 4));
  % z2 = find( (tcor >= 4) & (tcor < 10));
  % aa = nanmean(acor(z1));
  % anorm = nanmean(acor(z2));
  % if (anorm > 0)
  %    BI = (aa/anorm);
  % else
  %     BI = NaN;
  % end
  %*********
  
  %******* if you want to plot the result
  if (plotprogress)
    hf = figure(10);
    set(hf,'position',[1000 400 500 500]);
    plot(tcor,(acor ./ pcor),'k.-'); hold on;
    plot(tcor,ones(size(pcor)),'r-');
    plot(tcor,ncor,'b-');
    xlabel('Time (ms)');
    ylabel('Autocor');
    title(sprintf('BI=%5.2f',BI));
  end
  
return;


function smo = gauss_smooth(psth,Gsig)

    % Make the number of samples depending on the gaussian window size
    gaussian_filter_size = 4*Gsig-1; % if Gsig = 10, 19 samples total
                                     % 9 left & 9 right from the mean

    % Make smoothing kernel using gaussian filter
    for i = 1:gaussian_filter_size
        gauss  = exp(-(((i-(2*Gsig)).^2)./(2*Gsig^2)));
        gauss_filter(i,:) = gauss;
    end
    % Normalize the gaussian filter
    gauss_smooth = gauss_filter/sum(gauss_filter);
    psth_size    = length(psth);
    filter_size  = length(gauss_smooth);
    filter_cent = floor((filter_size+1)/2);

    for i=1:psth_size   % size_smooth

        % Always 0 for the initial value (only sum from product of two vectors)
        smo(i) = 0;
        nomo(i) = 0;

        % Apply filter to data
        for j = 1:filter_size
             diff = (j-filter_cent);   % this coordinate, goes - 3*sig up to 3*sig
             samp = (i+diff);
             if ( (samp >= 1) && (samp <= psth_size) )
                 smo(i) = smo(i) + (psth(samp) * gauss_smooth(j));
                 nomo(i) = nomo(i) + gauss_smooth(j);
             end       
        end
        %********
        if (nomo(i) > 0)
            smo(i) = smo(i) / nomo(i);
        end
        %***********
    end
    
return;

