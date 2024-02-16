function isp = SEC_makeSpline(prm,dataX)
    minx = min(dataX)/1.01;
    maxx = max(dataX)*1.01;

    bsp = spmak([minx minx prm(2) prm(2) maxx maxx],[prm(3), prm(4), prm(5)]);
    isp = fnint(bsp);
    isp = fnxtr(isp,2);
end