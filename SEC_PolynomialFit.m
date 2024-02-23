function [xbest,sigma] = SEC_PolynomialFit(polyDegree,dataX,dataY)

    % F = @(x) sum((10.^polyval(x(1:end),dataX)-dataY).^2);

    % xbest = polyfit(dataX,log10(dataY),polyDegree);
    xbest = polyfit(dataX,dataY,polyDegree);
    % %   Residuals
    % sigma = std(log10(dataY)-polyval(xbest,dataX));
    % x0 = [x0];
    % 
    % opt = optimset('Display','off');
    % xbest = fmincon(F,x0,[],[],[],[],[],[],[],opt);
end