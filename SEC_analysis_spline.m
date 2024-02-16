clear variables
close all
clc
%%  Data
%   In Excel
%   Column D - retention volume
%   Column E:N - 100nm, 80nm, 60nm, 50nm, 40nm, 30nm, 20nm, 15nm, 10nm, 5nm
chromatograms = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet', 'Chromatograms',...
    'Range','D5:N12004');

Rv = chromatograms(:,1);
%   Just based off what Lukas gave me
Rv = Rv(3301:7001);
Rv = round(Rv,5);
Extinction = chromatograms(3301:7001,2:end);


calib = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet','Calibration points',...
    'Range','C5:F13');

dataSizeAvg = calib(:,3);
dataSizeStdDev = calib(:,4);
dataVR = calib(:,2);

%% Fit to data

dataX = dataVR;
dataY = dataSizeAvg;
%   dataYvar = 1 if you consider all data points to have same uncertainty
%   otherwise provide variance of each point
dataYvar = dataSizeStdDev.^2;
% dataYvar = 1;

xbest = SEC_SplineFit(dataX,dataY,dataYvar);


tic
[samples,R,sampsThinned] = SEC_SplineStatisticalFit(dataX,dataY,dataYvar,...
    'burnin',0.5,...
    'nSamps',1000000,...
    'parallel',false,...
    'Nwalkers',100,...
    'stepSize',1.5);
toc

%%  Plot calibration curve fit


figure
%   Plot data + error bars
errorbar(dataVR,dataSizeAvg,dataSizeStdDev,'LineStyle','none','Marker','o','DisplayName','Data')
hold on
%   Plot fit from normal optimization
xplt = linspace(min(dataX),max(dataX));
isp = SEC_makeSpline(xbest,dataX);
yplt = xbest(1)+fnval(isp,xplt);
plot(xplt,yplt,'DisplayName','Spline Fit')

%   Plot statistical fit
indepSamps = sampsThinned(:,:);
nSimsStat = size(indepSamps,2);
plotSamps = indepSamps(:,randperm(nSimsStat));
statFit = zeros(nSimsStat,length(xplt));
for iii=1:nSimsStat
    isp = SEC_makeSpline(plotSamps(:,iii),dataX);
    statFit(iii,:) = plotSamps(1,iii) + fnval(isp,xplt);
end
Q = quantile(statFit,[0.005 0.025 0.5 0.975 0.995]);

plot(xplt,Q(3,:),...
    "Color",[0 0.4470 0.7410],...
    "LineStyle",'-',...
    'LineWidth',2,...
    'DisplayName','MCMC Median');

patch([xplt, fliplr(xplt)],...
    [Q(1,:), fliplr(Q(5,:))],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.1,...
    'DisplayName','MCMC 99%')

legend

hold off

%% Test
% function Y = testEval(x,prm,dataX)
%     minx = min(dataX)-1;
%     maxx = max(dataX)+1;
% 
%     bsp = spmak([minx minx prm(2) prm(2) maxx maxx],[prm(3), prm(4), prm(5)]);
%     isp = fnint(bsp);
%     isp = fnxtr(isp,2);
% 
%     Y = prm(1) + fnval(isp,x);
% end
% 
% function C = testCost(prm,dataX,dataY,dataYvar)
% 
%     simY = testEval(dataX,prm,dataX);
% 
%     C = sum( (dataY - simY).^2 ./dataYvar);
% 
% end

