clear variables
close all
clc
addpath("gwmcmc\")
rng('default')

%%  File to save workspace data to

svFileName = "polynomialCalibrationCurveFit.mat";
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

xbest = SEC_PolynomialFit(3,dataX,dataY,dataYvar);

tic
[~,R,sampsThinned] = SEC_PolynomialStatisticalFit(dataX,dataY,3,dataYvar,...
    'burnin',0.1,...
    'nSamps',10000000,...
    'parallel',false,...
    'Nwalkers',50,...
    'stepSize',1.5);
toc
%%
%   Do we "pass" the convergence tests
fprintf("___________________________________________________________\n")
fprintf("Rhat \t\t\t <1.01 is ideal\n")
fprintf("___________________________________________________________\n")
for iii=1:length(R)
    fprintf("Parameter %d \t\t %f\n",iii,R(iii))
end
fprintf("___________________________________________________________\n")

fprintf("# Independent Samples (min 100 recommended):\t %d\n",numel(sampsThinned/size(sampsThinned,1)))
fprintf("___________________________________________________________\n")

%% Plot calibration curve fit

%   Polynomial is fit to log of the data
plotX = linspace(min(dataX),max(dataX));
plotYlog10 = polyval(xbest,plotX);
plotY = 10.^plotYlog10;


%   Statistical fit
%   Use thinned samples to approximate independent draws from posterior
indepSamps = sampsThinned(:,:);
%   More samples are better but minimum of ~100 is ideal
% nSimsStat = size(indepSamps,2);
nSimsStat = 1000;
%   Randomly sort the samples array
plotSamps = indepSamps(:,randperm(nSimsStat));
%   Evaluate the polynomial for each sample
statFit = zeros(nSimsStat,length(plotX));
for iii=1:nSimsStat
    statFit(iii,:) = 10.^polyval(plotSamps(:,iii),plotX);
end
%   Quantiles of the statistical fit
Q = quantile(statFit,[0.005 0.025 0.5 0.975 0.995]);


figure
errorbar(dataVR,dataSizeAvg,dataSizeStdDev,'LineStyle','none','Marker','o','DisplayName','Data')
hold on
plot(plotX,plotY,'DisplayName','Normal Fit')

plot(plotX,Q(3,:),...
    "Color",[0 0.4470 0.7410],...
    "LineStyle",'-',...
    'LineWidth',2,...
    'DisplayName','MCMC Median');

patch([plotX, fliplr(plotX)],...
    [Q(1,:), fliplr(Q(5,:))],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.1,...
    'DisplayName','MCMC 99%')

legend
hold off


%%  Plot extinction-weighted PSD
nExtCurvs = size(Extinction,2);

%   Polynomial fit
xfitpoly = flip(10.^polyval(xbest,Rv));

qextpoly            = zeros(length(xfitpoly),nExtCurvs);
qextPolyErrorLo     = zeros(length(xfitpoly),nExtCurvs);
qextPolyErrorUp     = zeros(length(xfitpoly),nExtCurvs);
qextPolyMedian      = zeros(length(xfitpoly),nExtCurvs);

for iii=1:nExtCurvs
    qextpoly(:,iii) = xfitpoly .* flip(Extinction(:,iii));
    areapoly = trapz(xfitpoly,qextpoly(:,iii));
    qextpoly(:,iii) = qextpoly(:,iii) / areapoly;
    
    %   Bayesian fit median + 99% intervals
    qstatsplot = zeros(length(xfitpoly),nSimsStat);
    for jjj=1:nSimsStat
        xfitstats = flip(10.^polyval(plotSamps(:,jjj),Rv));
        qextstats = xfitstats .* flip(Extinction(:,iii));
        areastats = trapz(xfitstats,qextstats);
        qextstats = qextstats / areastats;
        qinterp = griddedInterpolant(xfitstats,qextstats,'makima','linear');
        qstatsplot(:,jjj) = qinterp(xfitpoly);
    end
    Qqext = quantile(qstatsplot,[0.005 0.5 0.995],2);

    qextPolyErrorLo(:,iii)      = Qqext(:,1);
    qextPolyErrorUp(:,iii)      = Qqext(:,3);
    qextPolyMedian(:,iii)       = Qqext(:,2);
    
    figure
    
    hold on
    plot(xfitpoly, qextpoly(:,iii) ,'DisplayName','Polyfit')
    
    plot(xfitpoly, qextPolyMedian(:,iii), 'DisplayName','Bayesian median')
    
    patch([xfitpoly; flip(xfitpoly)],...
        [qextPolyErrorLo(:,iii); flip(qextPolyErrorUp(:,iii))],...
        [0 0.4470 0.7410],...
        'FaceAlpha',0.1,...
        'DisplayName','Bayesian 99%')
    legend
    hold off
end

save(svFileName,'-mat')