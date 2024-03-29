function [samples,knot] = SEC_SplineStatisticalFit(dataX,dataY,dataYvar,options)

arguments
dataX {mustBeNonempty, mustBeFloat}
dataY {mustBeNonempty, mustBeFloat}
dataYvar {mustBeFloat} = 1
options.Nwalkers = 20;
options.startPerturb = 0.05
options.stepSize = 2.0
options.nSamps = 10000
options.thinChain = 1
options.burnin = 0.1
options.parallel = false
options.progressBar = true
end

%   Initial guess with normal optimization
[xfit,sigma2] = SEC_SplineFit(dataX,dataY);

%   Knot location is hard, so fix it to a single value now
knot = xfit(2);
nPrm = length(xfit)-1;

%   Likelihood based on sum of squares
%   exp(-0.5*sum of squares/variance) is form of normal distribution
%   log of this is just: -0.5*sum of squares/variance
logLikelihood = @(x) -0.5*SEC_SplineSumSquares([x(1);knot;x(2:end)],dataX,dataY,sigma2);

%   Use a flat prior distribution but constrain
% inBounds = @(prm) (prm(1)>0) * (prm(2)>=min(dataX)) * (prm(2)<=max(dataX)) * all(prm(3:end)<0);
inBounds = @(prm) (prm(1)>0) * all(prm(3:end)<0);
logPrior = @(x) -1/inBounds(x);

%   Log posterior is sum of prior and likelihood
% logPosterior = @(x) logPrior(x) + logLikelihood(x);
logPosterior = { logPrior logLikelihood};


%   Generate initial samples
xfitstat = [xfit(1) xfit(3:end)];
minit = (1+ options.startPerturb*unifrnd(-1,1,nPrm,options.Nwalkers)).*(xfitstat(:));
% lb = [0; min(dataX);-Inf;-Inf;-Inf];
% ub = [Inf; max(dataX);0;0;0]; 
lb = [0;-Inf;-Inf;-Inf];
ub = [Inf; 0;0;0]; 
for ccc=1:options.Nwalkers
    mcheck = minit(:,ccc);
    mcheck = min(mcheck,ub);
    mcheck = max(mcheck,lb);
    minit(:,ccc) = mcheck;
end

%   
samples = gwmcmc(minit, logPosterior, options.nSamps,...
    'StepSize',options.stepSize,...
    'ThinChain',options.thinChain,...
    'BurnIn',options.burnin,...
    'Parallel',options.parallel,...
    'ProgressBar',options.progressBar);

% R = psrf(samples);
% modelsFlat=models(:,:);
%
% modelsGR = permute(models,[3 1 2]);
% [~,~,ESS] = eacorr(samples);
% thin = (size(samples,2)*size(samples,3))/min(ESS);
% %   ESS ~100 at minimum as a rule of thumb
% 
% sampsThinned = samples(:,:,1:round(thin):end);
% % sampsThinnedFlat = sampsThinned(:,:);
% figure
% % Contours at 10%, 30%, 50%, 70%, 90%
% %
% ecornerplot(samples,'ks',true,'color',[0 0.4470 0.7410],'scatter',false,'grid',true)




end