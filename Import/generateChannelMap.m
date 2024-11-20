function indo = generateChannelMap(data,noise,start_indo)
% function indo = generateChannelMap(data,noise,start_indo)
    %inputs:
    %     data - the traces being fit
    %     noise - between 0 to 1, 0 always pick worst point to improve, 
    %               up to 1, choose them at random
    %     start_indo - use search start from previous index
    %
    %outputs:
    %     indo - best index order found
    %
    
%% plot the raw just to check that the channel map is correct

%inds = (1:30e3) + 55000e3+55000e3; % sample index
%inds = (30e3:1:60e3) + 55000e3+55000e3; % sample index another window
%inVolts = true;
%data = io.loadRaw(ops, inds, inVolts);
%data = bsxfun(@minus, data, mean(data)); % common average reference

Nchan = 64;
inda = [1:Nchan];
%********** here you would do an optimization of channel correlations
CMap = zeros(Nchan,Nchan);
for i = 1:Nchan
    disp(sprintf('Computing correlation chan %d',i));
    for j = i:Nchan
        [r,p] = corr(data(i,:)',data(j,:)');
        CMap(i,j) = r;
        CMap(j,i) = r;
    end
end
%************ search for superior ordering based on cost function *****
[UCost,Cost] = CorrCostFunction(CMap,inda);
OCost = Cost;
%******** search for better solution, print costs as you go *******
ICMap = zeros(Nchan,Nchan);
if (1)
  % indo = [1:33,58,39,34:35,40,36:38,41,50,43,42,44:47,49,51,53,48,55,52,57,54,56,59:64];
  if isempty(start_indo)
     indo = inda;
  else
     indo = start_indo;
  end
  for k = 1:40000
      [bindo,BCost,BUCost] = GetBestMove(CMap,indo,UCost,noise);
      disp(sprintf('Step %d: Cost %10.8f NCost %10.8f',k,OCost,BCost));
      if (BCost > Cost)
          Cost = BCost;
          indo = bindo;
          UCost = BUCost;
      end
      % input('check cost');
  end
  %***************
  idata(inda,:) = data(indo,:);  % remap the rows
  [IUCost,ICost] = CorrCostFunction(CMap,indo);
  %******* 
  for i = 1:Nchan
      for j = 1:Nchan
          ICMap(i,j) = CMap(indo(i),indo(j));
          ICMap(j,i) = ICMap(i,j);
      end
  end
  %*******
end
%***************

hf = figure; clf;
set(hf,'Position',[100 100 1200 900]);
%*******
subplot('Position',[0.1 0.1 0.25 0.8]);
channelOffsets = inda*200;
plot(bsxfun(@plus, data', channelOffsets)); hold on;   
for k = 1:length(indo)
    text(0.95*size(data,1),channelOffsets(k),sprintf('%d',indo(k)),'Fontsize',8);
end
title('Original Channel Map');
%******
subplot('Position',[0.4 0.6 0.25 0.30]);
imagesc(CMap,[0 1]); colorbar;
title(sprintf('Orig Cost:  %10.8f',OCost));
%****
subplot('Position',[0.4 0.1 0.25 0.30]);
imagesc(ICMap,[0 1]); colorbar;
title(sprintf('Rev Cost:  %10.8f',ICost));
%*******
subplot('Position',[0.7 0.1 0.25 0.8]);
plot(bsxfun(@plus, idata', channelOffsets)); hold on;
for k = 1:length(indo)
    text(0.95*size(data,1),channelOffsets(k),sprintf('%d',indo(k)),'Fontsize',10);
end
title('Revised Channel Map');
%*********************

return;


function [UnitCost,Cost,nindo] = CorrCostFunction(CMap,indo)

   N = size(CMap,1);
   UnitCost = zeros(1,N);
   for i = 1:N
       if (i == 1) || (i == 33)
           ucost = CMap(indo(i),indo(i+1));
       else
          if (i == 32) || (i == 64)
              ucost = CMap(indo(i),indo(i-1));
          else
              ucost = 0.5 * ( CMap(indo(i),indo(i-1)) + CMap(indo(i),indo(i+1)) );          
          end
       end
       UnitCost(i) = ucost;
   end
   Cost = sum(UnitCost)/N;  
return;
   
function [bindo,BCost,BUCost] = GetBestMove(CMap,indo,UnitCost,Noise)
   %******* find the best swap for worst costing unit ********
   N = size(CMap,1);
   UC = sort(UnitCost);
   Nscat = floor(Noise*N);
   if (Nscat < 1)
       Nscat = 1;
   end
   ii = randi(Nscat); % select from 3 worst options
   zworst = find( UnitCost == UC(ii));
   zw = zworst(1);
   %*********
   zindo = [1:(zw-1),(zw+1):N];
   ZN = length(zindo);
   %****** search for a better insertion point
   NuCost = zeros(1,N);
   BCost = 0;
   BUCost = UnitCost;
   if (zw <= 32)
       Zstart = 1;
       Zend = floor(N/2);
   else
       Zstart = floor(N/2)+1;
       Zend = N;
   end
   for i = 1:(ZN+1)
        if (i<Zstart) || (i>Zend)  % impose constrain shanks can't swap
            continue;
        end
        nindo = [zindo(1:(i-1)),zw,zindo(i:ZN)];
        [UCost,Cost] = CorrCostFunction(CMap,indo(nindo));
        if (Cost > BCost)
            BCost = Cost;
            bindo = indo(nindo);
            BUCost = UCost;
        end
   end
   %*************
return;
