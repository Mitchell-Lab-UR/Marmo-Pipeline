function PlotInfo_MoBeforeSaccade_PrefNPref_ROC(Info)
  
  %******** download variables stored in info
  fields = fieldnames(Info);
  for k = 1:size(fields,1)
      str = [fields{k} ' = Info.' fields{k} ';'];
      eval(str);
  end
  %*******************************
  
  %****** plot the results ****************
  H = figure;
  set(H,'Position',[100 100 300 900]);
  
  
  %********* PreSac - Plot AUC for Att and UnAtt
  subplot('position',[0.2 0.1 0.7 0.2])
  %****** subset of trials where no ori change occured
  xx = 1:(NPreBin-1);
  h2 = plot(TimeRoc(xx),ARoc(xx,1),'r-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(TimeRoc(xx),ARoc(xx,2),'r-');hold on;
  h2b = plot(TimeRoc(xx),ARoc(xx,3),'r-');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(TimeRoc(xx),URoc(xx,1),'b-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(TimeRoc(xx),URoc(xx,2),'b-');hold on;
  h2b = plot(TimeRoc(xx),URoc(xx,3),'b-');hold on;
  %************
  axis tight;
  V = axis;
  %axis([TimeRoc(1) TimeRoc(xx(end)) V(3) V(4)]);
  axis([-200 TimeRoc(xx(end)) V(3) V(4)]);
  plot([0,0],[V(3),V(4)],'k-');
  plot([TimeRoc(1),TimeRoc(xx(end))],[0.5,0.5],'k--');
  plot([PreWin(1),PreWin(1)],[V(3),V(4)],'k--');
  plot([PreWin(2),PreWin(2)],[V(3),V(4)],'k--');
  h = xlabel('Time (ms)');
  set(h,'FontSize',14);
  h = ylabel('AUC');
  set(h,'FontSize',14);
  h = title('Discrimination');
  set(h,'FontSize',14);
  
  %******** find max range *************
%   Vmax = max([max(AuuStimPref+(2*AsuStimPref)),max(AuuPreNPref+(2*AsuPrePref)),...
%               max(UuuStimPref+(2*UsuStimPref)),max(UuuPreNPref+(2*UsuPrePref))]);
  
  Vmax = 80; %300;
  %********* PreSac - attended, pref vs non-pref
  subplot('position',[0.2 0.7 0.7 0.2])
  %****** subset of trials where no ori change occured
  h2 = plot(-BefSac:AftSac,AuuPrePref,'g-');hold on;
  %set(h2,'Color',[1.0000 0.6000 0.7843])
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,AuuPrePref+(2*AsuPrePref),'g-');hold on;
  %set(h2b,'Color',[1.0000 0.6000 0.7843])
  h2b = plot(-BefSac:AftSac,AuuPrePref-(2*AsuPrePref),'g-');hold on;
  %set(h2b,'Color',[1.0000 0.6000 0.7843])
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,AuuPreNPref,'b-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,AuuPreNPref+(2*AsuPreNPref),'b-');hold on;
  h2b = plot(-BefSac:AftSac,AuuPreNPref-(2*AsuPreNPref),'b-');hold on;
  %************
  axis tight;
  V = axis;
  axis([-BefSac AftSac 0 Vmax]);
  plot([0,0],[0,Vmax],'k-');
  plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
  plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
  h = xlabel('Time (ms)');
  set(h,'FontSize',14);
  h = ylabel('Rate (hz)');
  set(h,'FontSize',14);
  h = title('Towards');
  set(h,'Color',[1,0,0]);
  set(h,'FontSize',14);
  
  %********* Stim - attended, pref vs non-pref
  subplot('position',[0.2 0.4 0.7 0.2])
  h2 = plot(-BefSac:AftSac,UuuPrePref,'g-');hold on;
  %set(h2,'Color',[ 0.5843 0.8157 0.9882])
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,UuuPrePref+(2*UsuPrePref),'g-');hold on;
  %set(h2b,'Color',[ 0.5843 0.8157 0.9882])
  h2b = plot(-BefSac:AftSac,UuuPrePref-(2*UsuPrePref),'g-');hold on; 
  %set(h2b,'Color',[ 0.5843 0.8157 0.9882])
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,UuuPreNPref,'b-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,UuuPreNPref+(2*UsuPreNPref),'b-');hold on;
  h2b = plot(-BefSac:AftSac,UuuPreNPref-(2*UsuPreNPref),'b-');hold on;
  %************
  axis tight;
  V = axis;
  axis([-BefSac AftSac 0 Vmax]);
  plot([0,0],[0,Vmax],'k-');
  plot([PreWin(1),PreWin(1)],[0,Vmax],'k--');
  plot([PreWin(2),PreWin(2)],[0,Vmax],'k--');
  h = xlabel('Time (ms)');
  set(h,'FontSize',14);
  h = ylabel('Rate (hz)');
  set(h,'FontSize',14);
  h = title('Away');
  h.Color = 'blue';
  set(h,'FontSize',14);
    
  return;
 