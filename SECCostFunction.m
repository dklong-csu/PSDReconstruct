function Cost = SECCostFunction(prm,splineOrder,nBasisFunctions,xMin,xMax,dataX,dataYMu, dataYSigma2)
nInteriorKnots = nBasisFunctions - splineOrder;

%   Make valid set of knots
knots = makeValidKnots(prm(1:nInteriorKnots),splineOrder,xMin,xMax);

%   Evaluate splines
sim = 0*dataX + prm(nInteriorKnots+1);
for iii=1:nBasisFunctions
    sim = sim + prm(nInteriorKnots+1+iii)*Ispline(dataX, splineOrder,iii,knots);
end

%   Cost function = (data - sim)^2/var
Cost = sum( (sim-dataYMu).^2 ./ dataYSigma2);
end