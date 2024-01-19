
k = 3;
left = 0;
right = 1;
interiorknots = [0.3, 0.5, 0.6];


t = [0.0,0.0,0.0, 0.3,0.5,0.6, 1.0,1.0,1.0];

x = linspace(0,1-1e-6);

M = zeros(6, length(x));
I = 0*M;
for iii=1:6
    M(iii,:) = Mspline(x,k,iii,t);
    I(iii,:) = Ispline(x,k,iii,t);
end
F = sum([1.2; 2.0; 1.2; 1.2; 3.0; 0.0] .* M,1);
F2 = sum([1.2; 2.0; 1.2; 1.2; 3.0; 0.0] .* I,1);

figure
hold on
for iii=1:6
    plot(x,M(iii,:));
    plot(x,I(iii,:))
    pause(1)
end
plot(x,F,'LineWidth',2)
plot(x,F2,'Linewidth',2)
hold off