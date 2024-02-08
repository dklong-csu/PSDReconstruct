%%  Fits SEC calibration curve and estimates the error

clc
clear variables
% close all
addpath("DRAM\")

rng('default')
%%  Todo
%   large retention volume -> size=0
%   continue slope for small retention volume?
%   Bayes factor
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

%%  Lukas method
fitvars = polyfit(dataVR,log10(dataSizeAvg),3);      
vr_plot = linspace(6,11,100);
y_fit = fitvars(1)*vr_plot.^3 + fitvars(2)*vr_plot.^2 + fitvars(3)*vr_plot + fitvars(4); 

%%  User settings
%   What degree polynomial?
k = 2;
%   How many interior points for the spline
n = 1;
%   How far beyond data should spline go?
factor = 1.01;
%   How many optimization iterations
nIters = 100;

%%   Implied information needed for spline
%   Minimum and Maximum retention volumes that need to be interpolated
VRmin = min(dataVR) / factor;
%   Add a small amount to VRMax because spline is not defined on the right
%   endpoint
VRmax = max(dataVR) * factor;

%   How many basis splines exist?
M = n+k;

%%  Fit calibration curve with no statistical information
%   Set the cost function
F = @(prm) SECCostFunction(prm,k,M,VRmin,VRmax,dataVR,dataSizeAvg,dataSizeStdDev.^2);
% F = @(prm) SECCostFunction(prm,k,M,VRmin,VRmax,dataVR,dataSizeAvg,1);

%   How many parameters are there
%   1 for each interiorKnot -- +n
%   1 for constant function -- +1
%   1 for coefficient of each basis -- +1
nPrm = n + 1 + M;

knots0 = linspace(VRmin,VRmax,n+2);
knots0 = knots0(2:end-1);
const0 = M;
coeff0 = -ones(1,M);
prm0 = [knots0, const0, coeff0];

lb = [(VRmin+1e-6)*(1+0*knots0), 0, (-Inf)*(1+0*coeff0)];
ub = [(VRmax-1e-6)*(1+0*knots0), Inf, 0*coeff0];
optOptions = optimset('Display','iter',...
    'MaxIter',nIters);
%   Inequality constraints
nConstr = n-1 + 1 + M + 1;
A = zeros(nConstr,length(prm0));
%   prm(1)<=prm(2)<=...<=prm(n)
%       knot(i+1) - knot(i) > 0.1
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
b(1:n-1) = -0.1 * ones(n-1,1);
bestPrmGlobal = simulannealbnd(F,prm0,lb,ub,optOptions);
bestPrm = fmincon(F,bestPrmGlobal,A,b,[],[],lb,ub,[],optOptions);

%%  Bayesian inversion
%   Evaluate log likelihood at parameters from fmincon


logLikeli = @(prm) -F(prm);
logPrior = @(prm) 1;

%   Negative Second deriv of log likelihood
clc
H = -computeHessian(logLikeli, bestPrm, 1e-3);
%   Laplace approximation: N(xbest, H^-1)
S = inv(H);

%   Initial points drawn from points near fmincon result
%       Rows are different parameters
%       Columns are different chains
n_chains = 4;
initS = 2.4^2/length(bestPrm) * S;
samp0 = mvnrnd(bestPrm,initS,n_chains)';
mcmc = DRAM(logLikeli, logPrior, samp0,...
    "initCov",initS,...
    "nChains",n_chains,...
    "nBurnedSamples",25,...
    "nAdaptSamples",0,...
    "nSamples",3000,...
    "updateEveryN",100,...
    "reprintHeaderEveryN",10,...
    "delayStages",0,...
    "covEpsilon",1,...
    "covarReduce",0.1);
mymap = parula;
vars = ["Knot1", "Knot2", "C0", "C1", "C2"];
t = plotPosterior(mcmc,vars,mymap);

%%  Plots

knots = makeValidKnots(bestPrm(1:n),k,VRmin,VRmax);
VRplot = linspace(VRmin,VRmax-1e-6);
diamSimPlot = 0*VRplot + bestPrm(n+1);
for iii=1:M
    diamSimPlot = diamSimPlot + bestPrm(n+iii+1)*Ispline(VRplot,k,iii,knots);
end

figure
errorbar(dataVR, dataSizeAvg, 2*dataSizeStdDev,...
    "vertical","x",...
    'DisplayName','Data',...
    'LineWidth',1.5,...
    'CapSize',10)


hold on
for iii=1:size(calib,1)
    %   FIXME - automate variable for the data
    mu = dataSizeAvg(iii);
    sigma = dataSizeStdDev(iii);
    nPlotPts = 1000;
    diamplot = linspace(mu-4*sigma,mu+4*sigma,nPlotPts);
    probs = normpdf(diamplot,mu, sigma);
    probs = probs/max(probs);
    plot(probs+dataVR(iii), diamplot,'HandleVisibility','off')
    
    patch([dataVR(iii)*ones(1,nPlotPts), fliplr(probs+dataVR(iii))],...
        [diamplot, fliplr(diamplot)],...
        'k',...
        'FaceAlpha',0.1,...
        'HandleVisibility','off')
end

plot(VRplot,diamSimPlot,'LineWidth',2,'DisplayName','Spline')
plot(vr_plot,10.^(y_fit),'LineStyle','--','LineWidth',2,'DisplayName','Cubic')
xlabel('Retention volume / ml')
ylabel('Particle core size / nm')
legend
% set(gca,'YScale','log')
hold off


%%
%--------------------------------------------------------------------------
%   Plot of polynomial basis functions
%--------------------------------------------------------------------------
figure
title('Fit of cubic to log10(data)')
ylabel('log10(diameter)')
xlabel('Retention volume')
scatter(dataVR, log10(dataSizeAvg),'DisplayName','Data')
hold on

y_fit3 = fitvars(1)*vr_plot.^3;
plot(vr_plot,y_fit3,'k--','DisplayName','Basis 1')

y_fit2 = fitvars(2)*vr_plot.^2;
plot(vr_plot,y_fit2,'k--','DisplayName','Basis 2')

y_fit1 = fitvars(3)*vr_plot.^1;
plot(vr_plot,y_fit1,'k--','DisplayName','Basis 3')

y_fit0 = fitvars(4)*vr_plot.^0;
plot(vr_plot,y_fit0,'k--','DisplayName','Basis 4')

plot(vr_plot,y_fit,'LineStyle','-','LineWidth',2,'DisplayName','Sum of Basis 1-4')
hold off



%%  Helper function
function H = computeHessian(F, x0, dx)
    %   Central differences
    ndims = length(x0);
    H = zeros(ndims,ndims);
    for iii=1:ndims
        for jjj=iii:ndims
            if iii==jjj
                h = dx * x0(iii);
                x1 = x0;
                x1(iii) = x1(iii)+h;
                x2 = x0;
                x2(iii) = x2(iii) - h;
                H(iii,iii) = ( F(x1) - 2*F(x0) + F(x2) )/(h^2);
            else
                hiii = dx * x0(iii);
                hjjj = dx * x0(jjj);
                
                x1 = x0;
                x1(iii) = x1(iii) + hiii;
                x1(jjj) = x1(jjj) + hjjj;
                f1 = F(x1);

                x2 = x0;
                x2(iii) = x2(iii) + hiii;
                x2(jjj) = x2(jjj) - hjjj;
                f2 = F(x2);

                x3 = x0;
                x3(iii) = x3(iii) - hiii;
                x3(jjj) = x3(jjj) + hjjj;
                f3 = F(x3);

                x4 = x0;
                x4(iii) = x4(iii) - hiii;
                x4(jjj) = x4(jjj) - hjjj;
                f4 = F(x4);

                hij = (f1 - f2 - f3 + f4)/(4*hiii*hjjj);
                H(iii,jjj) = hij;
                H(jjj,iii) = hij;
            end
        end
    end
end