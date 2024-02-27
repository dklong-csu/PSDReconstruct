function [xbest,sigma2] = SEC_SplineFit(dataX,dataY)
    F = @(p) SEC_SplineSumSquares(p,dataX,dataY,1);
    
    lb = [0, min(dataX),-Inf,-Inf,-Inf];
    ub = [Inf, max(dataX),0,0,0];   
    prm0 = [max(dataY),mean(dataX),-1,-1,-1];
    problem = createOptimProblem('fmincon',...
    'objective',F,...
    'x0',prm0,...
    'lb',lb,...
    'ub',ub);
    ms = MultiStart('Display','off','UseParallel',false);
    xbest = run(ms,problem,20);
    isp = SEC_makeSpline(xbest,dataX);
    simY = xbest(1) + fnval(isp,dataX);
    res = simY-dataY;
    sigma = std(res);
    sigma2 = sigma^2;
end
