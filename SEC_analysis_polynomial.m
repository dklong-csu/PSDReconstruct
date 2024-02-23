clear variables
close all
clc
addpath("gwmcmc\")
rng('default')

%%  Instructions
%       Make appropriate edits to the next two blocks "Program Controls"
%       and "Data"
%       The rest can be left alone unless specific edits are desired.

%%  Program Controls
polyDegree = 3;             %   what degree polynomial to fit to log10(diam)

%   For plotting error bars, how many samples should be used at maximum?
%   More is better, but takes longer
nSimsStatMax = 10000;

%   Generate samples?
makeSamples = true;
%   If false: assume you want to make plots of previous samples


%   MCMC general instructions
%       Try the recommended settings first based on polynomial degree
%       Check the summary statistics printed in the command window
%           -- If the Rhat table shows any value >1.01 --> increase nSamps
%              by a factor of 5 or 10
%           -- If the number of independent samples < 100 --> increase
%              nSamps approximately by a factor of 100/(independent
%              samples) or 500/(independent samples)

burn = 0.1;                 %   What percentage of samples to toss to 
                            %   "forget" the start of each chain?
                            %   (default=0.1)
nSamps=1e6;                 %   Total number of samples
                            %   degree=3 --> start with 1e6
                            %   degree=5 --> start with 1e7
Nwalkers=50;                %   Number of MCMC chains
                            %   degree=3 --> start with 50
stepSize=1.5;               %   >1, found 1.5 works fairly well for these problems

%   todo names for parameters

%   Generate text files to make plots directly in Latex?
%   This will create text files for each plot made with a column for the x
%   axis and a column for the y axis. So also useful if you want to make
%   plots in python or Origin or some other program
texPlots = true;

%   PSD plot diameters
xplot = linspace(0,150,200)';


%%  Data

%   Load the chromatograms that will be used to reconstruct PSDs
chromatograms = readmatrix('Calibration_core diameter.xlsx',...
    'Sheet', 'HPLC data',...
    'Range','H3:Q12003');

Rv = chromatograms(:,1);
Extinction = chromatograms(:,2:end);

%   Upper and lower bounds to filter out noise
%   Particularly important to set rvLo because the calibration curve
%   becomes large fast for small retention volumes, which is more of a
%   model limitation than an accurate result
rvLo = 5.8;
rvHi = 13;
%   Implements the filtering based on retention volume
idx = find(Rv > rvLo & Rv < rvHi);
Rv = Rv(idx);
Extinction = Extinction(idx,:);

%   Load the calibration data
%   Looking for a summary of
%       Data --> retention volume (x value) vs log10(diameter) (y value)
%       provide: 
%           mean value of log10(diameter) from TEM samples --> dataSizeAvg
%           sample variance of log10(diameter) from TEM    --> dataSizeVar
%           number of TEM samples                          -->dataSizeNTEM
%           peak retention volume corresponding to TEM samples --> dataVR
calib = readmatrix('Calibration_core diameter.xlsx',...
    'Sheet','HPLC data',...
    'Range','B14:E22');

dataSizeAvg = calib(:,2);
dataSizeVar = calib(:,3);
dataSizeNTEM = calib(:,4);
dataVR = calib(:,1);

%%  Files data are saved to 

%   Save workspace variables in this file
svFileName = strcat("degree",num2str(polyDegree),"CalibrationCurveFit.mat");

if ~makeSamples
    load(svFileName);
    makeSamples = false;
end

texPlotPrefix = strcat("degree",num2str(polyDegree),"Plots_");

%% Fit to data

if makeSamples
tic
samples = SEC_PolynomialStatisticalFit(dataVR,...
    dataSizeAvg,...
    dataSizeVar,...
    dataSizeNTEM,...
    polyDegree,...
    'burnin',burn,...
    'nSamps',nSamps,...
    'parallel',false,...
    'Nwalkers',Nwalkers,...
    'stepSize',stepSize);
toc
end

[sampsThinned, R, ESS] = SEC_analyze_samples(samples);
indepSamps = sampsThinned(:,:);
%   Plot statistics figure for samples
ecornerplot(samples,'ks',true,'color',[0 0.4470 0.7410],'scatter',false,'grid',true)
if texPlots
    writematrix(indepSamps',strcat(texPlotPrefix,"samples.txt"),'Delimiter','tab')
    fprintf("Independent Samples saved to: \t\t%s\n",strcat(texPlotPrefix,"samples.txt"))
end


%% Plot calibration curve fit

plotX = linspace(min(dataVR)/1.1,max(dataVR)*1.1);


%   Statistical fit
%   Use thinned samples to approximate independent draws from posterior
indepSamps = sampsThinned(:,:);
%   More samples are better but minimum of ~100 is ideal
nSimsStat = min(nSimsStatMax, size(indepSamps,2));
%   Randomly sort the samples array
plotSamps = indepSamps(:,randperm(nSimsStat));
%   Evaluate the polynomial for each sample
statFit = zeros(nSimsStat,length(plotX));
for iii=1:nSimsStat
    statFit(iii,:) = polyval(plotSamps(:,iii),plotX);
end
%   Quantiles of the statistical fit
Q = quantile(statFit,[0.005 0.025 0.5 0.975 0.995]);


figure
scatter(dataVR,dataSizeAvg,'Marker','o','DisplayName','Data')
% errorbar(dataVR,dataSizeAvg,dataSizeStdDev,'LineStyle','none','Marker','o','DisplayName','Data')
hold on
% plot(plotX,plotY,'DisplayName','Normal Fit')
% ylim([0 250])

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

if texPlots
    textCC = [plotX', Q(1,:)', Q(3,:)', Q(5,:)'];
    writematrix(textCC,strcat(texPlotPrefix,"CalibrationCurve.txt"),'Delimiter','tab')
    fprintf("Calibration curve saved to: \t\t%s\n",strcat(texPlotPrefix,"CalibrationCurve.txt"))
    fprintf("Format: \t\t\tRetention volume | Median log10(diameter) | 0.5%% log10(diameter) | 99.5%% log10(diameter)\n")
end


%%  Plot extinction-weighted PSD
nExtCurvs = size(Extinction,2);

%   Polynomial fit
%   Query relevant plotting diameters
% xplot = flip(10.^polyval(indepSamps(:,end),Rv));
% xplot = linspace(0,max(xplot)*1.2,200)';

qextpoly            = zeros(length(xplot),nExtCurvs);
qextPolyErrorLo     = zeros(length(xplot),nExtCurvs);
qextPolyErrorUp     = zeros(length(xplot),nExtCurvs);
qextPolyMedian      = zeros(length(xplot),nExtCurvs);

for iii=1:nExtCurvs
    %   Bayesian fit median + 99% intervals
    qstatsplot = zeros(length(xplot),nSimsStat);
    for jjj=1:nSimsStat
        xfitstats = flip(10.^polyval(plotSamps(:,jjj),Rv));
        qextstats = xfitstats .* flip(Extinction(:,iii));
        [xfitstats,III] = sort(xfitstats);
        qextstats = qextstats(III);
        areastats = trapz(xfitstats,qextstats);
        qextstats = qextstats / areastats;
        qinterp = griddedInterpolant([xfitstats(1)/1.01;xfitstats;xfitstats(end)*1.01],[0;qextstats;0],'linear','nearest');
        qstatsplot(:,jjj) = max(0,qinterp(xplot));
    end
    Qqext = quantile(qstatsplot,[0.005 0.5 0.995],2);

    qextPolyErrorLo(:,iii)      = Qqext(:,1);
    qextPolyErrorUp(:,iii)      = Qqext(:,3);
    qextPolyMedian(:,iii)       = Qqext(:,2);
    
    figure
    
    hold on
    
    plot(xplot, qextPolyMedian(:,iii), 'DisplayName','Bayesian median')
    
    patch([xplot; flip(xplot)],...
        [qextPolyErrorLo(:,iii); flip(qextPolyErrorUp(:,iii))],...
        [0 0.4470 0.7410],...
        'FaceAlpha',0.1,...
        'DisplayName','Bayesian 99%')
    legend
    hold off

    if texPlots
        psdTex = [xplot, Qqext(:,2), Qqext(:,1), Qqext(:,3)];
        writematrix(psdTex,strcat(texPlotPrefix,"PSD_",num2str(iii),".txt"),'Delimiter','tab')
        fprintf("PSD saved to: \t\t%s\n",strcat(texPlotPrefix,"PSD_",num2str(iii),".txt"))
        fprintf("Format: \t\t\tDiameter | Median qext | 0.5%% qext | 99.5%% qext\n")
    end
end

save(svFileName,'-mat')
fprintf("Workspace variables saved to: \t\t%s\n",svFileName)