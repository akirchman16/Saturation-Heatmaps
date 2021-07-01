clearvars;
close all;

% Uses SSA (Gillespie Algorithm) to simulate a simple reaction (A<-->B) and
% record the final population (set number of events) at the end of the
% simulation. A "heatmap" will then be generated based on the final
% population for various parameter sets.

k_f_Range = 0.1:0.1:10;   %range of values for k_f
k_r_Range = 0.1:0.1:10;   %range of values for k_r

Parameters = combvec(k_f_Range,k_r_Range);  %all parameter sets

% Memory Allocation
Final_A = zeros(1,length(Parameters));
Final_B = zeros(1,length(Parameters));

for Simulation = 1:length(Parameters)
    k_f = Parameters(1,Simulation); %value for k_f for this simulation
    k_r = Parameters(2,Simulation); %value for k_r for this simulation
    
    A = zeros(1,100); %memory allocation
    B = zeros(1,100); %memory allocation
    t = zeros(1,100);   %memory allocation
    
    A(1) = 100; %initial population for A
    B(1) = 0;   %initial population for B
    t(1) = 0;   %initial time
    
    Events = 0; %event counter for each simulation
    while Events <= 100 %amount of events occuring in each simulation
        Events = Events+1;  %advance event counter
        
        a = [k_f*A(Events), k_r*B(Events)];  %propensity functions for each reaction ([forward;backwards])
        a_0 = sum(a);   %total propensity function
        r = [rand, rand];    %random numbers for SSA
        dt = (1/a_0)*log(1/r(1));   %time to next reaction
        if a(1) > r(2)*a_0 %forward reaction occurs
            A(Events+1) = A(Events)-1;
            B(Events+1) = B(Events)+1;
        else   %reverse reaction occurs
            A(Events+1) = A(Events)+1;
            B(Events+1) = B(Events)-1;
        end
        
        t(Events+1) = t(Events)+dt; %advances time
    end
    Final_A(Simulation) = A(end);   %final populations
    Final_B(Simulation) = B(end);
end

[X,Y] = meshgrid(k_f_Range,k_r_Range);
Z_A = griddata(Parameters(1,:),Parameters(2,:),Final_A,X(:),Y(:));
Z_A = reshape(Z_A,size(X));
Z_B = griddata(Parameters(1,:),Parameters(2,:),Final_B,X(:),Y(:));
Z_B = reshape(Z_B,size(X));

figure(1);
imagesc(Parameters(1,:),Parameters(2,:),Z_A);
xlabel('k_f');
ylabel('k_r');
title('A Population - 100 Events');
colorbar('Ticks',0:10:100);