function xbest = SEC_SplineFit(dataX,dataY,dataYvar)
    F = @(p) SEC_SplineSumSquares(p,dataX,dataY,dataYvar);
    
    lb = [0, min(dataX),-Inf,-Inf,-Inf];
    ub = [Inf, max(dataX),0,0,0];   
    prm0 = [max(dataY),mean(dataX),-1,-1,-1];
    problem = createOptimProblem('fmincon',...
    'objective',F,...
    'x0',prm0,...
    'lb',lb,...
    'ub',ub);
    ms = MultiStart('Display','iter','UseParallel',false);
    xbest = run(ms,problem,20);
end
