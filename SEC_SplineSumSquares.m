function ssq = SEC_SplineSumSquares(prm, dataX, dataY, dataYvar)
    isp = SEC_makeSpline(prm,dataX);
    simY = prm(1) + fnval(isp,dataX);

    ssq = sum( (simY-dataY).^2 ./ dataYvar);

end