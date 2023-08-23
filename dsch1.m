function [tdsch,dsch1] = dsch1(r_ch,N)

n_i = round((length(N)-1)/4);   % node interval for cluster formation
dsch = zeros(40,4);

for j = 1:n_i
    dsch(j,1) = sqrt(((N(r_ch(1)).xp-N(j).xp)^2)+((N(r_ch(1)).yp-N(j).yp)^2));
end
tdsch(1,1) = sum(dsch(:,1));    % total distance between all the sensors to the cluster head in 1st cluster

for j = (n_i+1):(2*n_i)
    dsch(j,2) = sqrt(((N(r_ch(2)).xp-N(j).xp)^2)+((N(r_ch(2)).yp-N(j).yp)^2));
end
tdsch(1,2) = sum(dsch(:,2));    % total distance between all the sensors to the cluster head in 2nd cluster

for j = ((2*n_i)+1):(3*n_i)
    dsch(j,3) = sqrt(((N(r_ch(3)).xp-N(j).xp)^2)+((N(r_ch(3)).yp-N(j).yp)^2));
end
tdsch(1,3) = sum(dsch(:,3));    % total distance between all the sensors to the cluster head in 3rd cluster

for j = ((3*n_i)+1):(length(N)-1)
    dsch(j,4) = sqrt(((N(r_ch(4)).xp-N(j).xp)^2)+((N(r_ch(4)).yp-N(j).yp)^2));
end
tdsch(1,4) = sum(dsch(:,4));    % total distance between all the sensors to the cluster head in 4th cluster

dsch1 = zeros(40,1);
dsch1(1:10,1) = dsch(1:n_i,1);
dsch1(11:20,1) = dsch((n_i+1):(2*n_i),2);
dsch1(21:30,1) = dsch(((2*n_i)+1):(3*n_i),3);
dsch1(31:40,1) = dsch(((3*n_i)+1):(length(N)-1),4);

end