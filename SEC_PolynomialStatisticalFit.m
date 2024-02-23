function samples= SEC_PolynomialStatisticalFit(dataX,dataY,dataYvar,dataYN,polyDegree,options)

arguments
dataX {mustBeNonempty, mustBeFloat}
dataY {mustBeNonempty, mustBeFloat}
dataYvar {mustBeFloat} 
dataYN {mustBeNonempty, mustBeFloat}
polyDegree {mustBeInteger,mustBePositive} = 3
options.Nwalkers = 20;
options.startPerturb = 0.01
options.stepSize = 2.0
options.nSamps = 10000
options.thinChain = 1
options.burnin = 0.1
options.parallel = true
options.progressBar = true
end

%   Initial guess with normal optimization
xfit = SEC_PolynomialFit(polyDegree,dataX,dataY);
nPrm = length(xfit);

%   Likelihood based on sum of squares
%   exp(-0.5*sum of squares/variance) is form of normal distribution
%   log of this is just: -0.5*sum of squares/variance
logLikelihood = @(x) SEC_PolynomialLogStudentT(x,dataX,dataY,dataYvar,dataYN);
% logLikelihood = @(x) -0.5*sum((polyval(x(1:end),dataX)-log10(dataY)).^2./sigma.^2);

%   Use a flat prior distribution
% logPrior = @(x) -0.01*norm(x)^2;
% logPrior = @(x) log(1*all(abs(x)<1e12));
logPrior = @(x) 0;

%   Log posterior is sum of prior and likelihood
logPosterior = @(x) logPrior(x) + logLikelihood(x);


%   Generate initial samples
nInit = 0;
ntrys = 0;
minit = zeros(nPrm,options.Nwalkers);
while nInit < options.Nwalkers
    mtest = (1+ options.startPerturb*unifrnd(-1,1,nPrm,1)).*(xfit(:));
    testLL = logLikelihood(mtest);
    if isfinite(testLL)
        nInit = nInit + 1;
        minit(:,nInit) = mtest;
    end
    ntrys = ntrys+1;
end
% minit = (1+ options.startPerturb*unifrnd(-1,1,nPrm,options.Nwalkers)).*(xfit(:));

%   
samples = gwmcmc(minit, logPosterior, options.nSamps,...
    'StepSize',options.stepSize,...
    'ThinChain',options.thinChain,...
    'BurnIn',options.burnin,...
    'Parallel',options.parallel,...
    'ProgressBar',options.progressBar);

% R = psrf(samples);
% % modelsFlat=models(:,:);
% %
% % modelsGR = permute(models,[3 1 2]);
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