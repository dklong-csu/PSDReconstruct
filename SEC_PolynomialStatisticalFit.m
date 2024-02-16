function [samples,R,sampsThinned] = SEC_PolynomialStatisticalFit(dataX,dataY,polyDegree,dataYvar,options)

arguments
dataX {mustBeNonempty, mustBeFloat}
dataY {mustBeNonempty, mustBeFloat}
polyDegree {mustBeInteger,mustBePositive} = 3
dataYvar {mustBeFloat} = 1
options.Nwalkers = 20;
options.startPerturb = 0.05
options.stepSize = 2.0
options.nSamps = 10000
options.thinChain = 1
options.burnin = 0.1
options.parallel = true
options.progressBar = true
end

%   Initial guess with normal optimization
xfit = SEC_PolynomialFit(polyDegree,dataX,dataY,dataYvar);
nPrm = length(xfit);

%   Likelihood based on sum of squares
%   exp(-0.5*sum of squares/variance) is form of normal distribution
%   log of this is just: -0.5*sum of squares/variance
logLikelihood = @(x) -0.5*sum((10.^polyval(x,dataX)-dataY).^2./dataYvar);

%   Use a flat prior distribution
logPrior = @(x) 0;

%   Log posterior is sum of prior and likelihood
logPosterior = @(x) logPrior(x) + logLikelihood(x);


%   Generate initial samples
minit = (1+ options.startPerturb*unifrnd(-1,1,nPrm,options.Nwalkers)).*(xfit(:));

%   
samples = gwmcmc(minit, logPosterior, options.nSamps,...
    'StepSize',options.stepSize,...
    'ThinChain',options.thinChain,...
    'BurnIn',options.burnin,...
    'Parallel',options.parallel,...
    'ProgressBar',options.progressBar);

R = psrf(samples);
% modelsFlat=models(:,:);
%
% modelsGR = permute(models,[3 1 2]);
[~,~,ESS] = eacorr(samples);
thin = (size(samples,2)*size(samples,3))/min(ESS);
%   ESS ~100 at minimum as a rule of thumb

sampsThinned = samples(:,:,1:round(thin):end);
% sampsThinnedFlat = sampsThinned(:,:);
figure
% Contours at 10%, 30%, 50%, 70%, 90%
%
ecornerplot(samples,'ks',true,'color',[0 0.4470 0.7410],'scatter',false,'grid',true)




end