function knots = makeValidKnots(interiorKnots, splineOrder, xMin, xMax)
knots = [xMin*ones(1,splineOrder), interiorKnots, xMax*ones(1,splineOrder)];
end