function [ff] = FitnessFun(chr,num,den1,den2)
% chr - chromosome vale

dv = (bin2dec(chr))/10; % Decode value in 1 decimal point

ff = ((dv)*(num))/((den1)-(den2));

end