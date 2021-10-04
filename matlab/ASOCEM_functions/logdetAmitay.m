function logdet = logdetAmitay(A)
    [~,D] =  eig(A);
    D = diag(D);
    th = 10^-8;
    logdet = 0;
    for i=1:size(D,1)
        if D(i)>th
            logdet = logdet+log(D(i));
        end
    end
end
    