%%  Fits SEC calibration curve and estimates the error

clc
clear variables
close all

%%  Data
chromatograms = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet', 'Chromatograms',...
    'Range','D5:N12004');
%   In Excel
%   Column D - retention volume
%   Column E:N - 100nm, 80nm, 60nm, 50nm, 40nm, 30nm, 20nm, 15nm, 10nm, 5nm

calib_pts = readmatrix('Determination of calibration curve.xlsx',...
    'Sheet','Calibration points',...
    'Range','C5:F13');

%%  Make Plots

figure
errorbar(calib_pts(:,2), calib_pts(:,3), 2*calib_pts(:,4),...
    "vertical","o")


hold on
for iii=1:size(calib_pts,1)
    mu = calib_pts(iii,3);
    sigma = calib_pts(iii,4);
    diamplot = linspace(mu-4*sigma,mu+4*sigma);
    probs = normpdf(diamplot,mu, sigma);
    probs = probs/max(probs);
    plot(diamdistr+calib_pts(iii,2), diamplot)
    
    patch([calib_pts(iii,2)*ones(1,100), fliplr(diamdistr+calib_pts(iii,2))],...
        [diamplot, fliplr(diamplot)],...
        'k',...
        'FaceAlpha',0.1)
end
hold off