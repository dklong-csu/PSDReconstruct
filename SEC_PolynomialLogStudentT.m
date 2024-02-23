function LT = SEC_PolynomialLogStudentT(prm,dataRV,dataXmu,dataXsig,dataXn)

    % xsim = 10.^polyval(prm,dataRV);
    xsim = polyval(prm,dataRV);
    LT = 0;
    for iii=1:length(xsim)
        t = ( dataXmu(iii) - xsim(iii) )/sqrt(dataXsig(iii)/dataXn(iii));
        dof = dataXn(iii)-1;
        % LT = LT + log(tpdf(t,dof));
        %   Can exclude constant factors
        LT = LT - (dof+1)/2 * log(1 + t^2/dof);
    end
end