clc;
clear all;
close all;
tic

x_t = 100;  % considered length of x-axis
y_t = 100;  % considered length of y-axis
x_bs = 0.5*x_t; % x coordinate of base station
y_bs = 0.5*y_t; % y coordinate of base station
n_n = 40;  % number of nodes
n_i = round(n_n/4);    % node intervals to form 4 clusters in 4 quadrants
n_g = 20;   % number of generations
r = randperm(10,1);
r_ch = [r,r+10,r+20,r+30];
operating_nodes = n_n;
dead_nodes = 0;
temp_val = 0;
transmissions = 0;
flag1stdead = 0;
rnd = 1;  % Round of Operation

%%% Energy Values %%%
Eo = 2; % Initial Energy of a Node (in Joules)
% Energy required to run circuity (both for transmitter and receiver)
Eelec = 50*10^(-9); % units in Joules/bit
ETx = 50*10^(-9); % units in Joules/bit
ERx = 50*10^(-9); % units in Joules/bit
% Transmit Amplifier Types %
Eamp = 100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
% Data Aggregation Energy %
EDA = 5*10^(-9); % units in Joules/bit
% Size of data package %
k_bit = 4000; % units in bits
% random WSN creation

figure(1)
N(n_n+1).xp = x_bs;
N(n_n+1).yp = y_bs;
plot(N(n_n+1).xp,N(n_n+1).yp,'kx')
hold on
axis([0 100 0 100])

for i = 1:n_n
    
    N(i).cond = 1; % node is in operational if its value is 1 and non operational if its value is 0
    N(i).E = Eo; % alloting initial enery to nodes
    
    if (0<i) && (i<=n_i)
        N(i).xp = rand(1,1)*x_bs;    % sensor node x coordinate position
        N(i).yp = rand(1,1)*y_bs;    % sensor node y coordinate position
        N(i).type = 'SN';    % type of node
        if i == r_ch(1)
            N(i).type = 'CH';    % type of node
        end
    end
    
    if (n_i<i) && (i<=(2*n_i))
        N(i).xp = x_bs + rand(1,1)*x_bs;    % sensor node x coordinate position
        N(i).yp = rand(1,1)*y_bs;    % sensor node y coordinate position
        N(i).type = 'SN';    % type of node
        if i == r_ch(2)
            N(i).type = 'CH';    % type of node
        end
    end
    
    if ((2*n_i)<i) && (i<=(3*n_i))
        N(i).xp = rand(1,1)*x_bs;    % sensor node x coordinate position
        N(i).yp = y_bs + rand(1,1)*y_bs;    % sensor node y coordinate position
        N(i).type = 'SN';    % type of node
        if i == r_ch(3)
            N(i).type = 'CH';    % type of node
        end
    end
    
    if ((3*n_i)<i) && (i<=n_n)
        N(i).xp = x_bs + rand(1,1)*x_bs;    % sensor node x coordinate position
        N(i).yp = y_bs + rand(1,1)*y_bs;    % sensor node y coordinate position
        N(i).type = 'SN';    % type of node
        if i == r_ch(4)
            N(i).type = 'CH';    % type of node
        end
    end
    
    if N(i).type == 'CH'
        plot(N(i).xp,N(i).yp,'r*')
        hold on
    else
        plot(N(i).xp,N(i).yp,'bo')
        hold on
    end
    
    % encoding
    N(i).e = dec2bin(randperm(8,1)-1,3);    % node energy parameter
    N(i).d = dec2bin(randperm(8,1)-1,3);    % node distance parameter
    N(i).n = dec2bin(randperm(8,1)-1,3);    % node no. of neighbours parameter
    
    % chromosomes creation for initial population
    N(i).chr = strcat(N(i).e,N(i).d,N(i).n);
    
end
hold off

% total distance between all the sensors to the base station
dsbs = zeros(1,40);
for i = 1:n_n
    dsbs(i) = sqrt(((N(n_n+1).xp-N(i).xp)^2)+((N(n_n+1).yp-N(i).yp)^2));    % distance between sensor to the base station
end
tdsbs = sum(dsbs);


% generations process
while operating_nodes > 0
    
    energy = 0;
    dsch = zeros(1,40);
    
    figure(1)
    N(n_n+1).xp = x_bs;
    N(n_n+1).yp = y_bs;
    plot(N(n_n+1).xp,N(n_n+1).yp,'kx')
    hold on
    
    % total distance between all the sensors to the cluster head in respective clusters
    [tdsch,dsch] = dsch1(r_ch,N);
    
    % Energy Dissipation process
    for i = 1:n_n
        
        if (0<i) && (i<=n_i)
            % calculation of sensor node transmission energy
            if (N(i).E > 0) && (N(i).cond == 1) && (strcmp(N(i).type,'SN'))
                ETx = (Eelec * k_bit) + (Eamp * k_bit * (dsch(i))^2);
                N(i).E = N(i).E - ETx;
                energy = energy + ETx;
                ptfs(i) = ETx; % packets transmitted from source 
            elseif N(i).E <= 0  % if sensor node energy deplets
                N(i).cond = 0;
                N(i).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
            % calculation of cluster head receiver energy
            if (N(r_ch(1)).E > 0) && (N(r_ch(1)).cond == 1) && (strcmp(N(r_ch(1)).type,'CH'))
                ERx = (Eelec + EDA) * k_bit;
                N(r_ch(1)).E = N(r_ch(1)).E - ERx;
                energy = energy + ERx;
                ETx = ((Eelec + EDA) * k_bit) + (Eamp * k_bit * (dsbs(r_ch(1)))^2);
                N(r_ch(1)).E = N(r_ch(1)).E - ETx;
                energy = energy + ETx;
                prbs(i) = ETx; % packets received to base station
            elseif N(r_ch(1)).E <= 0    % if cluster head energy deplets
                N(r_ch(1)).cond = 0;
                N(r_ch(1)).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
        end
        
        if (n_i<i) && (i<=(2*n_i))
            % calculation of sensor node transmission energy
            if (N(i).E > 0) && (N(i).cond == 1) && (strcmp(N(i).type,'SN'))
                ETx = (Eelec * k_bit) + (Eamp * k_bit * (dsch(i))^2);
                N(i).E = N(i).E - ETx;
                energy = energy + ETx;
                ptfs(i) = ETx; % packets transmitted from source
            elseif N(i).E <= 0  % if sensor node energy deplets
                N(i).cond = 0;
                N(i).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
            % calculation of cluster head receiver energy and tranmission energy
            if (N(r_ch(2)).E > 0) && (N(r_ch(2)).cond == 1) && (strcmp(N(r_ch(2)).type,'CH'))
                ERx = (Eelec + EDA) * k_bit;
                N(r_ch(2)).E = N(r_ch(2)).E - ERx;
                energy = energy + ERx;
                ETx = ((Eelec + EDA) * k_bit) + (Eamp * k_bit * (dsbs(r_ch(2)))^2);
                N(r_ch(2)).E = N(r_ch(2)).E - ETx;
                energy = energy + ETx;
                prbs(i) = ETx; % packets received to base station
            elseif N(r_ch(2)).E <= 0    % if cluster head energy deplets
                N(r_ch(2)).cond = 0;
                N(r_ch(2)).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
        end

        if ((2*n_i)<i) && (i<=(3*n_i))
            % calculation of sensor node transmission energy
            if (N(i).E > 0) && (N(i).cond == 1) && (strcmp(N(i).type,'SN'))
                ETx = (Eelec * k_bit) + (Eamp * k_bit * (dsch(i))^2);
                N(i).E = N(i).E - ETx;
                energy = energy + ETx;
                ptfs(i) = ETx; % packets transmitted from source
            elseif N(i).E <= 0  % if sensor node energy deplets
                N(i).cond = 0;
                N(i).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
            % calculation of cluster head receiver energy
            if (N(r_ch(3)).E > 0) && (N(r_ch(3)).cond == 1) && (strcmp(N(r_ch(3)).type,'CH'))
                ERx = (Eelec + EDA) * k_bit;
                N(r_ch(3)).E = N(r_ch(3)).E - ERx;
                energy = energy + ERx;
                ETx = ((Eelec + EDA) * k_bit) + (Eamp * k_bit * (dsbs(r_ch(3)))^2);
                N(r_ch(3)).E = N(r_ch(3)).E - ETx;
                energy = energy + ETx;
                prbs(i) = ETx; % packets received to base station
            elseif N(r_ch(3)).E <= 0    % if cluster head energy deplets
                N(r_ch(3)).cond = 0;
                N(r_ch(3)).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
        end

        if ((3*n_i)<i) && (i<=n_n)
            % calculation of sensor node transmission energy
            if (N(i).E > 0) && (N(i).cond == 1) && (strcmp(N(i).type,'SN'))
                ETx = (Eelec * k_bit) + (Eamp * k_bit * (dsch(i))^2);
                N(i).E = N(i).E - ETx;
                energy = energy + ETx;
                ptfs(i) = ETx; % packets transmitted from source
            elseif N(i).E <= 0  % if sensor node energy deplets
                N(i).cond = 0;
                N(i).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
            % calculation of cluster head receiver energy
            if (N(r_ch(4)).E > 0) && (N(r_ch(4)).cond == 1) && (strcmp(N(r_ch(4)).type,'CH'))
                ERx = (Eelec + EDA) * k_bit;
                N(r_ch(4)).E = N(r_ch(4)).E - ERx;
                energy = energy + ERx;
                ETx = ((Eelec + EDA) * k_bit) + (Eamp * k_bit * (dsbs(r_ch(4)))^2);
                N(r_ch(4)).E = N(r_ch(4)).E - ETx;
                energy = energy + ETx;
                prbs(i) = ETx; % packets received to base station
            elseif N(r_ch(4)).E <= 0    % if cluster head energy deplets
                N(r_ch(4)).cond = 0;
                N(r_ch(4)).rop = rnd;
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
            end
        end
        
    end
    
    if (operating_nodes < n_n) && (temp_val == 0)
        toc
        temp_val = 1;
        flag1stdead = rnd;
    end
    
    op(rnd) = operating_nodes;
    
    if energy > 0
    nrg(rnd) = energy;
    end
        
    % fitness function calculation
    for i = 1:n_n
    
        if (0<i) && (i<=n_i)
            N(i).ffnum = sqrt(((N(r_ch(1)).xp-N(i).xp)^2)+((N(r_ch(1)).yp-N(i).yp)^2));
            N(i).ffden1 = tdsbs;
            N(i).ffden2 = tdsch(1,1);
            N(i).ff = FitnessFun(N(i).chr,N(i).ffnum,N(i).ffden1,N(i).ffden2);    % fitness function values for all sensors located in 1st cluster
        end

        if (n_i<i) && (i<=(2*n_i))
            N(i).ffnum = sqrt(((N(r_ch(2)).xp-N(i).xp)^2)+((N(r_ch(2)).yp-N(i).yp)^2));
            N(i).ffden1 = tdsbs;
            N(i).ffden2 = tdsch(1,2);
            N(i).ff = FitnessFun(N(i).chr,N(i).ffnum,N(i).ffden1,N(i).ffden2);   % fitness function values for all sensors located in 2nd cluster
        end

        if ((2*n_i)<i) && (i<=(3*n_i))
            N(i).ffnum = sqrt(((N(r_ch(3)).xp-N(i).xp)^2)+((N(r_ch(3)).yp-N(i).yp)^2));
            N(i).ffden1 = tdsbs;
            N(i).ffden2 = tdsch(1,3);
            N(i).ff = FitnessFun(N(i).chr,N(i).ffnum,N(i).ffden1,N(i).ffden2);   % fitness function values for all sensors located in 3rd cluster
        end

        if ((3*n_i)<i) && (i<=n_n)
            N(i).ffnum = sqrt(((N(r_ch(4)).xp-N(i).xp)^2)+((N(r_ch(4)).yp-N(i).yp)^2));
            N(i).ffden1 = tdsbs;
            N(i).ffden2 = tdsch(1,4);
            N(i).ff = FitnessFun(N(i).chr,N(i).ffnum,N(i).ffden1,N(i).ffden2);   % fitness function values for all sensors located in 4th cluster
        end

    end
    
    % probability density function
    pdf = ProDenFun(N);
    for i = 1:n_n
        N(i).pdf = pdf(i);
    end
    mpdf(1,1) = max(pdf(1:n_i));
    mpdf(1,2) = max(pdf((n_i+1):(2*n_i)));
    mpdf(1,3) = max(pdf(((2*n_i)+1):(3*n_i)));
    mpdf(1,4) = max(pdf(((3*n_i)+1):n_n));
    
    % assigning new cluster heads for the next generation
    for i = 1:n_n
        
        if (0<i) && (i<=n_i)
            if N(i).pdf == mpdf(1,1)
                N(i).type = 'CH';    % assigning this node as cluster head
                r_ch(1) = i;
            else
                N(i).type = 'SN';
            end
        end

        if (n_i<i) && (i<=(2*n_i))
            if N(i).pdf == mpdf(1,2)
                N(i).type = 'CH';    % assigning this node as cluster head
                r_ch(2) = i;
            else
                N(i).type = 'SN';
            end
        end

        if ((2*n_i)<i) && (i<=(3*n_i))
            if N(i).pdf == mpdf(1,3)
                N(i).type = 'CH';    % assigning this node as cluster head
                r_ch(3) = i;
            else
                N(i).type = 'SN';
            end
        end

        if ((3*n_i)<i) && (i<=n_n)
            if N(i).pdf == mpdf(1,4)
                N(i).type = 'CH';    % assigning this node as cluster head
                r_ch(4) = i;
            else
                N(i).type = 'SN';
            end
        end
        
        if N(i).type == 'CH'
            plot(N(i).xp,N(i).yp,'r*')
            hold on
        else
            plot(N(i).xp,N(i).yp,'bo')
            hold on
        end
        
    end
    
    % Crossover operation to generate new population for next generation
    i = 1;
    A = 1:n_n;
    while (i <= (n_n/2))

        B = A(randperm(length(A),1));
        A(A==B) = [];
        C = A(randperm(length(A),1));
        A(A==C) = [];

        D = strcat(N(B).chr(1:5),N(C).chr(6:9));
        E = strcat(N(C).chr(1:5),N(B).chr(6:9));

        N(B).chr = D;
        N(C).chr = E;

        i = i+1;

    end
    
    tprbs(rnd) = sum(prbs);
    tptfs(rnd) = sum(ptfs);
    
    pdr1(rnd) = tprbs(rnd)/tptfs(rnd);

    hold off
    rnd = rnd+1;
    
end
toc

figure(2)
plot(op,'-r','Linewidth',2)
pdr = flip(pdr1);
title ({'GA based clustering algorithm'; 'Operating Nodes per Round';})
xlabel 'Rounds';
ylabel 'Operational Nodes';
figure(3)
plot(pdr)
title ({'GA based clustering algorithm'; 'Packet Delivery Ratio per Round';})
xlabel 'Rounds';
ylabel 'Packet Delivery Ratio';
figure(4)
plot(nrg)
title ({'GA based clustering algorithm'; 'Energy consumption per Round';})
xlabel 'Rounds';
ylabel 'Energy';