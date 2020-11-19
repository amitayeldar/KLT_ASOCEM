function logdet = logdetAmitay(A)
    th = 10^-8;
    [~,D] =  eig(A);
    D = diag(D);
    logdet = 0;
    for i=1:size(D,1)
        if D(i)>th
            logdet = logdet+log(D(i));
        end
    end
end
    