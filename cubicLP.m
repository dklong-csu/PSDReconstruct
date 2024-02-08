function [lp, glp] = cubicLP(prm,dataX,dataY,dataStd)
    nprm = length(prm);
    nx = length(dataX);
    lp = 0;
    glp = zeros(nprm,1);
    for iii=1:nx
        y = polyval(prm,dataX(iii));
        d = y-dataY(iii);
        v = dataStd(iii)^2;
        
        lp = lp - d^2/v;
        glp = glp - 2/v * d * [dataX(iii)^3;dataX(iii)^2;dataX(iii);1];
    end

end