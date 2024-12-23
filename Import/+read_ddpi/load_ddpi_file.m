function vpx  = load_ddpi_file( filename )
% function  vpx = load_ddpi_file( filename )
% - reads ddpi file out for raw data, median filters and smooths
%   to form a file up-sampled to 1000 hz
% - also detects trial start and end strobe times with taglets
%***
%  input:  filename
%***
%  outputs:   vpx struct
%             vpx.raw - [Nx4] - fields: time, x, y, pupil
%             vpx.smo - [Tx4] - same fields, but upsampled 1000hz
%                               via linear interp, and smoothing with
%                               a Gaussian of sigma 5ms
%             returns [] if file problem
%
%***** it reads the VPX file strobes and reformats them
%***** to be just like ephys strobes, so we can treat them the same
%             vpx.tstrobes - timestamps of strobes
%             vpx.strobes - strobe values
%****** HOWEVER, it reads from a DDPI file format and converts to
%****** a format consistent with previous VPX struct
%*******************************************

MEDFILT = 0;   % leave raw data from DDPI (no median filtering)
% DDPI is more clean and Jake may want the raw
%*******************
fprintf('Reading DDPI file %s\n',filename);
vpx = struct;
vpx_strobes = struct;

%***** read in the raw eye data first
Output = read_ddpi.ddpiReadFile(filename);
% 1) signalType;
% 2) time;
% 3) p1x;
% 4) p1y;
% 5) p1r;
% 6) p1I;
% 7) p4x;
% 8) p4y;
% 9) p4r;
% 10) p4I;
% 11) p4score;
% 12) tag;
% 13) message;

%%
tt = Output(2,:)';
tt = tt / 1000;  % covert ms into secs to be consistent with VPX

% P1
p1x = Output(3,:)';
p1y = Output(4,:)';
p1r = Output(5,:)';
p1I = Output(6,:)';
p4x = Output(7,:)';
p4y = Output(8,:)';
n = numel(p1x);

% detect good tracking
figure(1); clf
subplot(1,2,1)
plot(p1x, p4x, '.');
axis tight
hold on
xlabel('P1')
ylabel('P4')
title('Horizontal Position')

mdlrx = fitlm(p1x,p4x,'RobustOpts','on');
mdlry = fitlm(p1y,p4y,'RobustOpts','on');
xcoef = mdlrx.Coefficients.Variables;
ycoef = mdlry.Coefficients.Variables;

xd = xlim;
plot(xd, xd*xcoef(2,1) + xcoef(1,1), 'r')

subplot(1,2,2)
plot(p1y, p4y, '.');
axis tight
hold on
xlabel('P1')
ylabel('P4')
title('Vertical Position')
xd = xlim;
plot(xd, xd*ycoef(2,1) + ycoef(1,1), 'r')

% use residuals to find outliers
p4xhat = p1x*xcoef(2,1) + xcoef(1,1);
residx = p4x - p4xhat;

p4yhat = p1y*ycoef(2,1) + ycoef(1,1);
residy = p4y - p4yhat;

outliers = false(n,1);
for i = 1:5
    z = hypot(zscore(residx(~outliers)),zscore(residy(~outliers)));
    outliers(~outliers) = z > 4.5;
%     sum(outliers)
end

figure(2); clf
subplot(1,2,1)
plot(p1x, p4x, '.')
hold on
plot(p1x(~outliers), p4x(~outliers), '.')
xlabel('P1')
ylabel('P4')
title('X position')

subplot(1,2,2)
plot(p1y, p4y, '.')
hold on
plot(p1y(~outliers), p4y(~outliers), '.')
xlabel('P1')
ylabel('P4')
title('Y position')
legend({'outliers', 'good'})

figure(3); clf
gazex = p4x - p1x;

plot(gazex); hold on
plot(find(~outliers), gazex(~outliers), '.')
xlabel('Sample #')
ylabel('Position')

figure(4); clf
plot(p1x(~outliers), p4x(~outliers) - p4xhat(~outliers), '.'); hold on
plot(xlim, [0 0], 'k--')

% %%
% X = [p1x(~outliers) p1x(~outliers).^2];
% Y = [p1y(~outliers) p1y(~outliers).^2];
% mdlrx = fitlm(X,p4x(~outliers),'RobustOpts','on');
% mdlry = fitlm(Y,p4y(~outliers),'RobustOpts','on');
% xcoef = mdlrx.Coefficients.Variables;
% ycoef = mdlry.Coefficients.Variables;
% 
% p4xhat = X*xcoef(2:end,1) + xcoef(1,1);
% 
% 
% %%
% % figure(1); clf
% % plot(p1x(~outliers), p4xhat, '.')
% % 
% % v = [p1x(~outliers), p4xhat];
% u = [p1x(~outliers) p4x(~outliers)];
% % 
% % xx = linspace(min(p1x(~outliers)), max(p1x(~outliers)), 1e3)';
% % yy = [xx(:) xx(:).^2]*xcoef(2:end,1) + xcoef(1);
% % plot(xx, yy)
% % 
% % d = hypot(u(:,1)-xx', u(:,2)-yy');
% % [~, id] = min(d, [], 2);
% % 
% % p4x2 = yy(id);
% % p1x2 = xx(id);
% % %%
% figure(1); clf
% 
% n = size(u,1);
% ns = 10;
% starts = [1:ceil(n/ns):n n];
% cmap = hsv(ns);
% for i = 1:ns
%     plot(u(starts(i):starts(i+1),1), u(starts(i):starts(i+1),2), '.', 'Color', cmap(i,:)); hold on
% end
% plot(p1x, p4x, '.'); hold on
% plot(p1x2, p4x2, '.'); hold on

% % plot(u(:,2) - u(:,1))
%%
ex = p4x - p1x;
ey = p4y - p1y;

ex(outliers) = nan;
ey(outliers) = nan;




%% --- find borked P1 times
% speed of intensity or radius changes
p1m = hypot([0; diff(p1I)], [0; diff(p1r)]);
% no intensity change for multiple frames is indicative of a borked p1
p1bork = imgaussfilt(p1m, 2)==0;

% running median intensity
T = 540*10;
m = medfilt1(p1I, T, 'truncate');
% 50% of the running intensity
thresh = .5*m;

p1bork = p1bork | p1I < thresh;
p1bork = imgaussfilt(double(p1bork), 2)>0;

blink = p1bork;

ep = 1e3*ones(size(ex));
ep(blink) = 0;

dat = [tt(:) ex(:) ey(:) ep(:)];
vpx.raw = dat;


%***** RUN a median filter on the data that is +/- 2 samples ********
disp('Median filtering plus and minus 2 samples');
mdat = dat;  % smooth traces remaining
if (MEDFILT > 0)
    fprintf('Median Filtering Eye Data (+/- %d samples)\n',MEDFILT);
    for k = (MEDFILT+1):(size(dat,1)-MEDFILT-1)
        mdat(k,2) = median( dat((k-MEDFILT):(k+MEDFILT),2) );
        mdat(k,3) = median( dat((k-MEDFILT):(k+MEDFILT),3) );
        mdat(k,4) = median( dat((k-MEDFILT):(k+MEDFILT),4) );
        %**********
        if (mod(k,10000) == 0)
            fprintf('Med filt DDPI rawdata, count %d\n',k);
        end
        %*********
    end
else
    fprintf('Running SGOLAYFILT on DDPI data ...\n');
    mdat(:,2) = sgolayfilt(mdat(:,2),1,3);
    mdat(:,3) = sgolayfilt(mdat(:,3),1,3);
end
%******** for now, only median filtered for smoothed
vpx.smo = mdat;
disp('Completed DDPI smoothing and upsampling');

%***************** now return to the file and build up strobe history
[vpx.tstrobes, vpx.strobes] = read_ddpi.read_ddpi_strobes( filename );
%************************************************************
start = find(Output(1,:)==1, 1, 'first');
stop = find(Output(1,:)==1, 1, 'last');
vpx.smo = vpx.smo(start:stop,:);
vpx.raw = vpx.raw(start:stop,:);


