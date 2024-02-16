function xbest = SEC_PolynomialFit(polyDegree,dataX,dataY,dataYvar)

    F = @(x) sum((10.^polyval(x,dataX)-dataY).^2./dataYvar);

    nDims = polyDegree + 1;

    x0 = polyfit(dataX,log10(dataY),polyDegree);

    xbest = fmincon(F,x0);

end