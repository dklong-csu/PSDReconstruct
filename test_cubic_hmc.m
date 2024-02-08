%%  Fits SEC calibration curve and estimates the error

clc
clear variables
close all

rng('default')
%%  Data
chromatograms = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet', 'Chromatograms',...
    'Range','D5:N12004');
%   In Excel
%   Column D - retention volume
%   Column E:N - 100nm, 80nm, 60nm, 50nm, 40nm, 30nm, 20nm, 15nm, 10nm, 5nm

calib = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet','Calibration points',...
    'Range','C5:F13');

dataSizeAvg = calib(:,3);
dataSizeStdDev = calib(:,4);
dataVR = calib(:,2);

%%  
for iii=1:length(dataVR)
    fitY(iii) = log10( dataSizeAvg(iii)^2 / sqrt(dataSizeAvg(iii)^2 + dataSizeStdDev(iii)^2));
    fitStd(iii) = sqrt(log10(1 + dataSizeStdDev(iii)^2/dataSizeAvg(iii)^2));
    fitX(iii) = dataVR(iii);
end

x0 = [1,1,1,1];
lpost = @(p) cubicLP(p,fitX,fitY,fitStd);

smp = hmcSampler(lpost,x0);

xhat = estimateMAP(smp,'VerbosityLevel',1)
tsmp = tuneSampler(smp);

nChains=4;
chains = cell(nChains,1);
eps = 0.1;
nsamps = 1000;
nburn = 500;

for iii=1:nChains
    pert = unifrnd(-1,1,length(xhat),1);
    s0 = xhat + eps * pert.*xhat;
    chains{iii} = drawSamples(tsmp,'Burnin',nburn,'NumSamples',nsamps,...
        'Start',s0,...
        'VerbosityLevel',1,...
        'NumPrint',100);
end
chainstats = diagnostics(tsmp,chains)
%%
figure
errorbar(dataVR,dataSizeAvg,2*dataSizeStdDev,"linestyle",'none')
hold on
yvals = polyval(xhat,fitX);
yvals = 10.^(yvals);
plot(fitX,yvals,'LineWidth',2)
for iii=1:50
    c = randi(nChains,1);
    s = randi(nsamps,1);
    a = chains{c};
    p = a(s,:);
    yvals = polyval(p, fitX);
    yvals = 10.^(yvals);
    plot(fitX,yvals,'k--')

end
hold off