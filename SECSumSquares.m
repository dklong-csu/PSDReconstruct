function Cost = SECSumSquares(prm,splineOrder,nBasisFunctions,xMin,xMax,dataX,dataYMu)
nInteriorKnots = nBasisFunctions - splineOrder;

%   Make valid set of knots
knots = makeValidKnots(prm(1:nInteriorKnots),splineOrder,xMin,xMax);

%   Evaluate splines
sim = 0*dataX + prm(nInteriorKnots+1);
for iii=1:nBasisFunctions
    sim = sim + prm(nInteriorKnots+1+iii)*Ispline(dataX, splineOrder,iii,knots);
end

Cost = (sim-dataYMu).^2 ;
end