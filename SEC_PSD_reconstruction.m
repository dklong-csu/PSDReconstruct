%%  Fits SEC calibration curve and estimates the error

clc
clear variables
close all

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

%%  Standard way to fit calibration curve
%   What degree polynomial?
k = 2;
%   How many interior points for the spline
n = 2;

%   Implied information
%   Minimum and Maximum retention volumes that need to be interpolated
% VRmin = min([chromatograms(:,1); calib(:,2)]);
VRmin = min(calib(:,2));
%   Add a small amount to VRMax because spline is not defined on the right
%   endpoint
% VRmax = max([chromatograms(:,1); calib(:,2)]) + 1e-6;
VRmax = max(calib(:,2)) + 1;

%   How many basis splines exist?
M = n+k;

%   Set the cost function
% F = @(prm) SECCostFunction(prm,k,M,VRmin,VRmax,calib(:,2),calib(:,3),calib(:,4).^2);
F = @(prm) SECCostFunction(prm,k,M,VRmin,VRmax,calib(:,2),calib(:,3),1);

prm0 = [8,12,90,-10,-10,-10,-10,-10];
lb = [VRmin+1e-6,VRmin+1e-6,0,-Inf,-Inf,-Inf,-Inf,-Inf];
ub = [VRmax-1e-6,VRmax-1e-6,Inf,0,0,0,0,0];
optOptions = optimset('Display','iter',...
    'MaxFunEvals',2000);
%   Inequality constraints
nConstr = n-1 + 1 + M + 1;
A = zeros(nConstr,length(prm0));
%   prm(1)<=prm(2)<=...<=prm(n)
for iii=1:n-1
    A(iii,iii:iii+1) = [1 -1];
end
%   prm(n+1)>0
A(n,n+1) = -1;
%   prm(n+2:end)<=0
A(n+1:length(prm0)-1, n+2:length(prm0)) = eye(length(n+2:length(prm0)));
%   sum(prm(n+1:end))>0
A(end,n+1:end) = -ones(1,length(prm0)-n);

b = zeros(size(A,1),1);
bestPrmGlobal = simulannealbnd(F,prm0,lb,ub,optOptions);
bestPrm = fmincon(F,bestPrmGlobal,A,b,[],[],lb,ub,[],optOptions);

%%  Plots
% bestPrm = [10,12,18,64.7416,-15.8706,3.5989,-15.4297,-26.1301,-5.1722,13.6363];

knots = makeValidKnots(bestPrm(1:n),k,VRmin,VRmax);
VRplot = linspace(VRmin,VRmax-1e-6);
diamSimPlot = 0*VRplot + bestPrm(n+1);
for iii=1:M
    diamSimPlot = diamSimPlot + bestPrm(n+iii+1)*Ispline(VRplot,k,1,knots);
end
figure
errorbar(calib(:,2), calib(:,3), 2*calib(:,4),...
    "vertical","o")


hold on
for iii=1:size(calib,1)
    mu = calib(iii,3);
    sigma = calib(iii,4);
    diamplot = linspace(mu-4*sigma,mu+4*sigma);
    probs = normpdf(diamplot,mu, sigma);
    probs = probs/max(probs);
    plot(probs+calib(iii,2), diamplot)
    
    patch([calib(iii,2)*ones(1,100), fliplr(probs+calib(iii,2))],...
        [diamplot, fliplr(diamplot)],...
        'k',...
        'FaceAlpha',0.1)
end

plot(VRplot,diamSimPlot)
hold off