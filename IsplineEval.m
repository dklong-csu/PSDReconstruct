function Y = IsplineEval(x,prm,splineOrder,nBasisFunctions,xMin,xMax)

nInteriorKnots = nBasisFunctions - splineOrder;
knots = makeValidKnots(prm(1:nInteriorKnots),splineOrder,xMin,xMax);
Y = 0*x + prm(nInteriorKnots+1);
for iii=1:nBasisFunctions
    Y = Y + prm(nInteriorKnots+1+iii)*Ispline(x, splineOrder,iii,knots);
end

end