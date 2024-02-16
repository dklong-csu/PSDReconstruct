%%  Fits SEC calibration curve and estimates the error

clc
clear variables
% close all
addpath("gwmcmc\")

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
nIters = 500;

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
[bestPrm,fvalmin,~,~,~,df,ddf] = fmincon(F,bestPrmGlobal,A,b,[],[],lb,ub,[],optOptions);

%%  Bayesian inversion
%   Evaluate log likelihood at parameters from fmincon
%   Knot location is super sensitive and difficult to infer so let us just
%   set it as constant
iKnot = bestPrm(1:n);
% 
logLikeli = @(prm) -SECCostFunction([iKnot(:);prm(:)],k,M,VRmin,VRmax,dataVR,dataSizeAvg,dataSizeStdDev.^2);
%   Prior distribution matching linear constrains in fmincon. This will
%   result in "1" if constraints are met (i.e. log prior=-1) and "0" if
%   constraints are not met (i.e. log prior = -Inf)
logPrior = @(prm) -1 ./ ( all(A*([iKnot(:);prm(:)]) <= b) .* all([iKnot(:);prm(:)] <= ub') .* all([iKnot(:);prm(:)] >= lb') );

logPfun = @(prm) logLikeli(prm) + logPrior(prm);
Nwalker=100;
minit = (1+ 0.02*unifrnd(-1,1,nPrm-n,Nwalker)).*(bestPrm(n+1:end)');

tic
rng('default')
nSamps = 100000;
models = gwmcmc(minit, logPfun, nSamps,...
    'StepSize',1.5,...
    'ThinChain',1,...
    'BurnIn',0.1,...
    'Parallel',true);
toc
R = psrf(models)
modelsFlat=models(:,:);
%
% modelsGR = permute(models,[3 1 2]);
[~,~,ESS] = eacorr(models)
thin = (size(models,2)*size(models,3))/min(ESS)
%   ESS ~100 at minimum as a rule of thumb

sampsThinned = models(:,:,1:round(thin):end);
sampsThinnedFlat = sampsThinned(:,:);
figure
% Contours at 10%, 30%, 50%, 70%, 90%
%
ecornerplot(models,'ks',true,'color',[0 0.4470 0.7410],'scatter',false,'grid',true)

%%
VRplot = linspace(VRmin,VRmax-1e-6)';
modelfcn = @(x,prm) IsplineEval(x,prm,k,M,VRmin,VRmax);




postfcn = zeros(size(sampsThinnedFlat,2),length(VRplot));
tic
parfor iii=1:size(sampsThinnedFlat,2)
    postfcn(iii,:) = modelfcn(VRplot,[iKnot;sampsThinnedFlat(:,iii)]);
end
toc
Q = quantile(postfcn,[0.005 0.025 0.5 0.975 0.995]);
%%
% for iii=1:size(models,2)
%     plot(VRplot, modelfcn(VRplot, models(:,iii)),'k')
% end
figure
hold on
plot(VRplot,Q(3,:),...
    "Color",[0 0.4470 0.7410],...
    "LineStyle",'-',...
    'LineWidth',2,...
    'DisplayName','MCMC Median');

patch([VRplot', fliplr(VRplot')],...
    [Q(1,:), fliplr(Q(5,:))],...
    [0 0.4470 0.7410],...
    'FaceAlpha',0.1,...
    'DisplayName','MCMC 99%')




%  Plots

knots = makeValidKnots(bestPrm(1:n),k,VRmin,VRmax);
VRplot = linspace(VRmin,VRmax-1e-6);
diamSimPlot = 0*VRplot + bestPrm(n+1);
for iii=1:M
    diamSimPlot = diamSimPlot + bestPrm(n+iii+1)*Ispline(VRplot,k,iii,knots);
end

% figure
errorbar(dataVR, dataSizeAvg, 2*dataSizeStdDev,...
    "vertical","x",...
    'DisplayName','Data',...
    'LineWidth',1.5,...
    'CapSize',10)


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






