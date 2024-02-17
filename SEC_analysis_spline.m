clear variables
close all
clc
rng('default')
addpath("gwmcmc\")
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
[samples,R,sampsThinned,knot] = SEC_SplineStatisticalFit(dataX,dataY,dataYvar,...
    'burnin',0.1,...
    'nSamps',1000000,...
    'parallel',false,...
    'Nwalkers',100,...
    'stepSize',1.5);
toc
R
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
    prm = [plotSamps(1,iii);knot;plotSamps(2:end,iii)];
    isp = SEC_makeSpline(prm,dataX);
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


%%  Plot extinction-weighted PSD

%   Bayesian fit median + 99% intervals
qinterpstats = cell(nSimsStat,1);
xmin = Inf;
xmax = 0;
for iii=1:nSimsStat
    prm = [plotSamps(1,iii);knot;plotSamps(2:end,iii)];
    isp = SEC_makeSpline(prm,dataX);
    %   Spline has trouble with extrapolation so do it myself
    lb = isp.breaks(2);
    ub = isp.breaks(end-1);

    Rvmid = Rv(and(Rv>=lb,Rv<=ub));
    xvalmid = plotSamps(1,iii) + fnval(isp,Rvmid);

    %   Find secant line between second two knots and use that as slope for
    %   extrapolation
    Rvlo = Rv(Rv<lb);
    m = ( fnval(isp,isp.breaks(3)) - fnval(isp,isp.breaks(2)) )/ ( isp.breaks(3) - isp.breaks(2) );
    x1 = lb;
    y1 = plotSamps(1,iii) + fnval(isp,x1);
    %   y-y1 = m(x-x1)
    xvallo = m*(Rvlo-x1)+y1;

    %   Secant line from second to last pair of knots to extrapolate
    Rvhi = Rv(Rv>ub);
    m = ( fnval(isp,isp.breaks(end-2)) - fnval(isp,isp.breaks(end-1)) )/( isp.breaks(end-2) - isp.breaks(end-1) );
    x1 = ub;
    y1 = plotSamps(1,iii) + fnval(isp, x1);
    %   y-y1 = m(x-x1)
    xvalhi = m*(Rvhi-x1)+y1;

    %   Combine interp and extrap
    xvalsstat = [xvallo; xvalmid; xvalhi]; % need to flip
    xvalsstat = flip(xvalsstat);
    %   Remove negative sizes
    goodidx = xvalsstat>0;
    goodxval = xvalsstat(goodidx);
    extin = flip(Extinction(:,1));
    goodExtinction = max(0,extin(goodidx));

    %   PSD
    qextstat = goodxval .* goodExtinction;
    area = trapz(goodxval,qextstat);
    qextstat = qextstat/area;
    qinterpstats{iii} = griddedInterpolant(goodxval,qextstat,'makima','nearest');
    xmin = min(xmin,min(goodxval));
    xmax = max(xmax,max(goodxval));
end

xvals = linspace(xmin,xmax,1000)';
qstatsplot = zeros(length(xvals),nSimsStat);
for iii=1:nSimsStat
    qinterp = qinterpstats{iii};
    qstatsplot(:,iii) = qinterp(xvals);
end
Qqext = quantile(qstatsplot,[0.005 0.025 0.5 0.975 0.995],2);

figure

hold on
% plot(xfitpoly, qextpoly,'DisplayName','Polyfit')

plot(xvals, Qqext(:,3), 'DisplayName','Bayesian median')

patch([xvals; flip(xvals)],...
    [Qqext(:,1); flip(Qqext(:,5))],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.1,...
    'DisplayName','Bayesian 99%')
legend
hold off