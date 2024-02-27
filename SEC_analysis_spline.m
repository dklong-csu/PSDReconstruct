clear variables
close all
clc
rng('default')
addpath("gwmcmc\")

%%  Program Controls

%   Save workspace variables in this file
svFileName = "splineCalibrationCurveFitMany.mat";

%   For plotting error bars, how many samples should be used at maximum?
%   More is better, but takes longer
nSimsStatMax = 1000;

%   Generate samples?
makeSamples = false;
%   If false: assume you want to make plots of previous samples
%   What is the filename to load for that?
if ~makeSamples
    load("splineCalibrationCurveFitMany.mat");
    makeSamples = false;
end

%   todo controls here for mcmc settings

%   todo names for parameters

%   Generate text files to make plots directly in Latex?
%   This will create text files for each plot made with a column for the x
%   axis and a column for the y axis. So also useful if you want to make
%   plots in python or Origin or some other program
texPlots = true;
texPlotPrefix = "splinePlots_";

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

% xbest = SEC_SplineFit(dataX,dataY,dataYvar);

if makeSamples
    tic
    [samples,knot] = SEC_SplineStatisticalFit(dataX,dataY,dataYvar,...
        'burnin',0.1,...
        'nSamps',10000000,...
        'parallel',false,...
        'Nwalkers',100,...
        'stepSize',1.5);
    toc
end

[sampsThinned, R, ESS] = SEC_analyze_samples(samples);
ecornerplot(samples,'ks',true,'color',[0 0.4470 0.7410],'scatter',false,'grid',true)
if texPlots
    indepSamps = sampsThinned(:,:);
    writematrix(indepSamps',strcat(texPlotPrefix,"samples.txt"),'Delimiter','tab')
end

%%  Plot calibration curve fit


figure
%   Plot data + error bars
errorbar(dataVR,dataSizeAvg,dataSizeStdDev,'LineStyle','none','Marker','o','DisplayName','Data')
hold on
%   Plot fit from normal optimization
xplt = linspace(min(dataX)/1.1,max(dataX)*1.1);
% isp = SEC_makeSpline(xbest,dataX);
% yplt = xbest(1)+fnval(isp,xplt);
% plot(xplt,yplt,'DisplayName','Spline Fit')

%   Plot statistical fit
indepSamps = sampsThinned(:,:);
nSimsStat = min(nSimsStatMax, size(indepSamps,2));
plotSamps = indepSamps(:,randperm(nSimsStat));
statFit = zeros(nSimsStat,length(xplt));
for iii=1:nSimsStat
    prm = [plotSamps(1,iii);knot;plotSamps(2:end,iii)];
    isp = SEC_makeSpline(prm,dataX);
    %   Spline has trouble extrapolating so do it myself
    %   small RV/large diam extrap
    m1 = ( fnval(isp,isp.breaks(3)) - fnval(isp,isp.breaks(2)) )/ ( isp.breaks(3) - isp.breaks(2) );
    smallRV = min(xplt);
    xsmallRV = m1*(smallRV - isp.breaks(2)) + plotSamps(1,iii) + fnval(isp,isp.breaks(2));
    %   large RV/small diam extrap
    m2 = ( fnval(isp,isp.breaks(end-2)) - fnval(isp,isp.breaks(end-1)) )/( isp.breaks(end-2) - isp.breaks(end-1) );
    largeRV = max(xplt);
    xlargeRV = m2*(largeRV - isp.breaks(end-1)) + plotSamps(1,iii) + fnval(isp,isp.breaks(end-1));
    %   interpolation object to eval easier
    rvInterp = linspace(min(dataX),max(dataX),1000);
    xInterp = plotSamps(1,iii) + fnval(isp,rvInterp);
    ccExtrap = griddedInterpolant([smallRV,rvInterp,largeRV],[xsmallRV,xInterp,xlargeRV],'makima','linear');
    % statFit(iii,:) = plotSamps(1,iii) + fnval(isp,xplt);
    statFit(iii,:) = ccExtrap(xplt);
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

if texPlots
    textCC = [xplt', Q(1,:)', Q(3,:)', Q(5,:)'];
    writematrix(textCC,strcat(texPlotPrefix,"CalibrationCurve.txt"),'Delimiter','tab')
end


%%  Plot extinction-weighted PSD
nExtCurvs = size(Extinction,2);

for jjj=1:nExtCurvs
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
        extin = flip(Extinction(:,jjj));
        goodExtinction = max(0,extin(goodidx));
    
        %   PSD
        qextstat = goodxval .* goodExtinction;
        area = trapz(goodxval,qextstat);
        qextstat = qextstat/area;
        qinterpstats{iii} = griddedInterpolant(goodxval,qextstat,'makima','nearest');
        xmin = min(xmin,min(goodxval));
        xmax = max(xmax,max(goodxval));
    end
    
    xvals = linspace(xmin,xmax,200)';
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

    if texPlots
        psdTex = [xvals, Qqext(:,3), Qqext(:,1), Qqext(:,5)];
        writematrix(psdTex,strcat(texPlotPrefix,"PSD_",num2str(jjj),".txt"),'Delimiter','tab')
    end
end

save(svFileName,'-mat')