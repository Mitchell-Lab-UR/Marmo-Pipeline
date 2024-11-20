function PlotInfo_MoBeforeSaccade(Info)
  
  %******** download variables stored in info
  fields = fieldnames(Info);
  for k = 1:size(fields,1)
      str = [fields{k} ' = Info.' fields{k} ';'];
      eval(str);
  end
  %*******************************
  
  %****** plot the results ****************
  H = figure;
  set(H,'Position',[100 100 300 300]);
   
  %*****************************************
  subplot('position',[0.2 0.2 0.7 0.7]);
  %*** plot fits over raw data points
  zYJit = max(max(Atunepre.mu),max(Utunepre.mu)) * YJit;
  PlotWithVonMisesDots(OriVals,TrA,zPreOri,PreSpk{NPreBin},Atunepre,'r',XJit,zYJit);
  PlotWithVonMisesDots(OriVals,TrU,zPreOri,PreSpk{NPreBin},Utunepre,'b',XJit,zYJit);
  %********** loop around
  if (1)
      OriVals = [OriVals;(OriVals(1)+360)];
      Atunepre.mu = [Atunepre.mu Atunepre.mu(1)];
      Atunepre.sem = [Atunepre.sem Atunepre.sem(1)];
      Utunepre.mu = [Utunepre.mu Utunepre.mu(1)];
      Utunepre.sem = [Utunepre.sem Utunepre.sem(1)];
  end
  %****** plot vonMises fit
  h2 = plot(OriVals,Atunepre.mu + (2 * Atunepre.sem),['r-']);
  h2 = plot(OriVals,Atunepre.mu - (2 * Atunepre.sem),['r-']);
  h2 = plot(OriVals,Utunepre.mu + (2 * Utunepre.sem),['b-']);
  h2 = plot(OriVals,Utunepre.mu - (2 * Utunepre.sem),['b-']);
  %**********
  h2 = plot(OriVals,Atunepre.mu,['r-']);
  set(h2,'Linewidth',2);
  h2 = plot(OriVals,Utunepre.mu,['b-']);
  set(h2,'Linewidth',2);
  %***********
  maxo = max(max(Atunepre.mu+(2*Atunepre.sem)),max(Utunepre.mu+(2*Utunepre.sem))) * 1.2;
  mino = min(min(Atunepre.mu+(2*Atunepre.sem)),min(Utunepre.mu+(2*Utunepre.sem))) * 0.0;
  axis([0 360 mino maxo]);
  h = xlabel('Direction (degs)');
  set(h,'FontSize',14);
  h = ylabel('Spike Counts');
  set(h,'FontSize',14);
  h = title('Tuning');
  set(h,'FontSize',14);
  %************  
   
  return;
  
  
 function PlotWithVonMisesDots(OriVals,TrA,zPreOri,StimSpk,Atune,colo,XJit,YJit)
  
   %******
   OriNum = size(OriVals,1);
   for i = 1:OriNum
        %****** plot raw counts
        % first plot cases where saccade to RF
        pzz = find( (zPreOri(TrA) == OriVals(i)));
        xjit = (rand(size(pzz))-0.5)*XJit;
        yjit = (rand(size(pzz))-0.5)*YJit;     
        %********
        H = plot(OriVals(i)*ones(size(pzz))+xjit,StimSpk(TrA(pzz))+yjit,[colo,'.']); hold on;
        set(H,'Markersize',5);
        if (colo == 'r')
            set(H,'Color','w');
            %set(H,'Color',[1,0.6,0.6]);
        end
        if (colo == 'b')
            set(H,'Color','w');
            %set(H,'Color',[0.6,0.6,1]);
        end 
        
  end
  
return;