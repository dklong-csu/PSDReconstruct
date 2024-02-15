function R = psrf(samples)

%   samples - nPrm x nChains x nSamps
[D,C,N] = size(samples);
split = zeros(D,2*C,floor(N/2)); % fixme
split(:,1:C,:) = samples(:,:,1:floor(N/2));
split(:,C+1:2*C,:) = samples(:,:,floor(N/2)+1:2*floor(N/2));

R = zeros(1,D);
for ddd=1:D
    d = squeeze(split(ddd,:,:));
    meanall = mean(d,"all");
    meanchain = mean(d,2);
    B = (N/2)/(2*C-1) * sum( (meanchain-meanall).^2 );
    W = mean(var(d,0,2));

    n = floor(N/2);
    varp = (n-1)/n * W + B/n;
    R(ddd) = sqrt(varp/W);
end

end