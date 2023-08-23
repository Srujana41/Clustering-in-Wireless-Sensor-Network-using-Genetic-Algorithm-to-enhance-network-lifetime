function [pdf] = ProDenFun(N)

    n_i = round((length(N)-1)/4);   % node interval for cluster formation
    
    sff = 0;
    for j = 1:n_i
        sff = sff + N(j).ff;
    end
    sffc(1,1) = sff;    % sum of fitness function values in 1st cluster

    sff = 0;
    for j = (n_i+1):(2*n_i)
        sff = sff + N(j).ff;
    end
    sffc(1,2) = sff;    % sum of fitness function values in 2nd cluster

    sff = 0;
    for j = ((2*n_i)+1):(3*n_i)
        sff = sff + N(j).ff;
    end
    sffc(1,3) = sff;    % sum of fitness function values in 3rd cluster

    sff = 0;
    for j = ((3*n_i)+1):(length(N)-1)
        sff = sff + N(j).ff;
    end
    sffc(1,4) = sff;    % sum of fitness function values in 4th cluster
    
    pdf = zeros(1,(length(N)-1));
    for i = 1:(length(N)-1)

        if (0<i) && (i<=n_i)
            pdf(i) = N(i).ff/sffc(1,1);   % probability density function values for all sensors located in 1st cluster
        end

        if (n_i<i) && (i<=(2*n_i))
            pdf(i) = N(i).ff/sffc(1,2);   % probability density function values for all sensors located in 2nd cluster
        end

        if ((2*n_i)<i) && (i<=(3*n_i))
            pdf(i) = N(i).ff/sffc(1,3);   % probability density function values for all sensors located in 3rd cluster
        end

        if ((3*n_i)<i) && (i<=(length(N)-1))
            pdf(i) = N(i).ff/sffc(1,4);   % probability density function values for all sensors located in 4th cluster
        end

    end

end