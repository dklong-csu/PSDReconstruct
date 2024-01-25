%%
% close all
% clear vars
% clc
%%  Settings
interiorKnots = bestPrm(1:n);
splineDegree = k;
leftEndPoint = 6.1370/1.01;
rightEndPoint = 10.2835*1.01;
constFcn = bestPrm(n+1);
coeffs = bestPrm(n+2:end);

%%  Form spline
knots = makeValidKnots(interiorKnots, splineDegree, leftEndPoint, rightEndPoint);
nBasisFcns = length(interiorKnots) + splineDegree;


x = linspace(leftEndPoint,rightEndPoint-1e-6,1000);

M = zeros(nBasisFcns, length(x));
I = 0*M;
for iii=1:nBasisFcns
    M(iii,:) = Mspline(x,splineDegree,iii,knots);
    I(iii,:) = Ispline(x,splineDegree,iii,knots);
end
F = sum(coeffs .* M,1);
F2 = constFcn + sum(coeffs .* I,1);

%%  Plot spline
figure
hold on
for iii=1:nBasisFcns
    % plot(x,M(iii,:));
    plot(x,I(iii,:))
    pause(1)
end
% plot(x,F,'LineWidth',2)
plot(x,F2,'Linewidth',2)
hold off