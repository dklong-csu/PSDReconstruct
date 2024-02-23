function [samplesThinned, R, ESS] = SEC_analyze_samples(samples)

R = psrf(samples);
%   Get median and 99% range
Q = quantile(samples(:,:),[0.005 0.5 0.995],2);

%   Do we "pass" the convergence tests
sep_str = "____________________________________________________________________________________________\n";
fprintf("\nModel:\t log10(diameter)\t=\t")
for iii=1:length(R)
    fprintf("C%d * RV^%d",iii,length(R)-iii)
    if iii<length(R)
        fprintf("\t+\t")
    end
end
fprintf("\n")
fprintf(sep_str)
fprintf("Rhat \t\t\t <1.01 is ideal\t\t\t(0.5%% \t\t--\t\t 50%% \t\t--\t\t 99.5%%)\n")
fprintf(sep_str)
for iii=1:length(R)
    fprintf("C%d \t\t %f \t\t %10.6e \t--\t %10.6e \t--\t %10.6e\n",iii,R(iii),Q(iii,1),Q(iii,2),Q(iii,3))
end
fprintf(sep_str)



[~,~,ESS] = eacorr(samples);
thin = (size(samples,2)*size(samples,3))/min(ESS);
%   ESS ~100 at minimum as a rule of thumb

samplesThinned = samples(:,:,1:round(thin):end);

fprintf("# Total Samples                            :\t %d\n",numel(samples)/size(samples,1))
fprintf("# Independent Samples (min 100 recommended):\t %d\n",numel(samplesThinned)/size(samplesThinned,1))
fprintf(sep_str)

end