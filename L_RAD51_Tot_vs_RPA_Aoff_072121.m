clearvars
close all;

% Generates a heatmap figure for variations of the two parameters
% L_RAD51_Total and k_off_RPA_A. It allows for plotting the growth 
% profile of each simulation but also the heatmaps. Heatmaps generated 
% include RAD51 saturation, RPA saturation, and time to equilibrium.

tic

N = 5000;   %length of ssDNA lattice
% RAD51 Parameters
RAD51 = 51; %value that will be stored on lattice to represent bound RAD51
n_RAD51 = 3;    %length of RAD51 protein

L_RAD51_Total_Values = linspace(1,10,2);  %total concentration of RAD51 in solution
Percent_M_RAD51 = 0.5;    %percentage of RAD51 solution which is monomers
w_RAD51 = 1;    %cooperativity parameter for RAD51
k_on_RAD51 = 1;     %kinetic rate constant for RAD51 binding to ssDNA
k_off_RAD51 = 1;    %kinetic rate constant for RAD51 dissociating from ssDNA

% RPA Parameters
RPA_A = 1;  %value to represent A piece of RPA on lattice
RPA_D = 3;  %value to represent D piece of RPA on lattice
n_A = 10;   %length of A component of RPA
n_D = 10;   %length of D component of RPA

L_RPA = 2;  %concentration of RPA in solution
w_RPA = 1;  %cooperativity parameter of RPA (for macroscopic binding)
k_on_RPA_A = 25; %kinetic rate constant for RPA-A binding to ssDNA
k_on_RPA_D = 15;  %kinetic rate constant for RPA-D binding to ssDNA
k_off_RPA_A_Values = linspace(1,10,2); %kinetic rate constant for RPA-A dissociating from ssDNA
k_off_RPA_D = 1; %kinetic rate constant for RPA-D dissociating from ssDNA

n_RPA = sum([n_A,n_D]);   %calculates total length of RPA molecule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Heatmap_Parameters = combvec(L_RAD51_Total_Values,k_off_RPA_A_Values);   %combination of possible sets of parameters that will generate heatmap (for this instance: L_RAD51_Total and k_off_RPA_A)

% Memory Allocation
RPA_Avg_Saturation = zeros(1,numel(Heatmap_Parameters)/2);
RAD51_Avg_Saturation = zeros(1,numel(Heatmap_Parameters)/2);
Total_Avg_Saturation = zeros(1,numel(Heatmap_Parameters)/2);
t_Equilibrium = zeros(1,numel(Heatmap_Parameters)/2);
Simulation_Times = zeros(1,numel(Heatmap_Parameters)/2);
% Allocation of cell arrays for data that will produce growth profiles
time_Data = cell(numel(Heatmap_Parameters)/2,1);
FracCover_RAD51_Data = cell(numel(Heatmap_Parameters)/2,1);
FracCover_RPA_A_Data = cell(numel(Heatmap_Parameters)/2,1);
FracCover_RPA_D_Data = cell(numel(Heatmap_Parameters)/2,1);
FracCover_RPA_Data = cell(numel(Heatmap_Parameters)/2,1);
FracCover_Total_Data = cell(numel(Heatmap_Parameters)/2,1);

P1 = Heatmap_Parameters(1,:);  %Temp. variable (for L_RAD51_Total) to avoid high, over-necessary communication in parfor of Parameters
P2 = Heatmap_Parameters(2,:);  %Temp. variable (for k_off_RPA_A) to avoid high, over-necessary communication in parfor of Parameters

w = waitbar(0,['Running ', num2str(numel(Heatmap_Parameters)/2), ' Simulations...']);  %generates a progress bar for the simulations
q_Count = parallel.pool.DataQueue;  %data queue to count how many simulations are completed
afterEach(q_Count,@parforWaitbar);  %runs parforWaitbar function at conclusion of each iteration (updates progress bar)
parforWaitbar(w, numel(Heatmap_Parameters)/2); %defines inputs for parforWaitbar function

parfor Simulations = 1:(numel(Heatmap_Parameters)/2)
    tic
    % Referencing each parameter value for each simulation
    L_RAD51_Total = P1(Simulations); %total concentration of RAD51
    k_off_RPA_A = P2(Simulations);    %kinetic rate constant for unbinding of RPA-A
    
    L_RAD51_M = L_RAD51_Total*Percent_M_RAD51;  %calculates concentration of RAD51 monomers
    L_RAD51_D = L_RAD51_Total-L_RAD51_M;    %calculates concentration of RAD51 dimers

    % Memory Allocation
    DNA = zeros(2,N);   %array to represent DNA
                        %top row is used to store locations of hinged open RPA-D
                        %bottom row actually represents the DNA itself
    RAD51_Mon_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Monomers are bound
    RAD51_Dim_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Dimers are bound
    RPA_A_BoundAtSpot = zeros(1,N); %array to record where RPA-A is actively bound
    RPA_D_BoundAtSpot = zeros(1,N); %array to record where RPA-D is actively bound
    RPA_D_HingedOpen = zeros(1,N);  %array to record where RPA-D is microscopically dissociated from lattice
    LocationHistory = zeros(14,1);  %Matrix used to store locations of all events. Same order as Full_Propensity

    % Initial Values
    t = 0;   %initial time is zero
    FracCover_RAD51 = 0; %initially no RAD51 is on the lattice
    FracCover_RPA_A = 0; %initially empty lattice
    FracCover_RPA_D = 0; %initially empty lattice
    FracCover_RPA = 0;   %initially no RPA on the lattice
    FracCover_Total = 0; %initially the lattice is empty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Event_Count = 0;    %counts how many events happen within the simulation
    Equilibrium_RAD51 = 0;  %test of whether RAD51 saturation is at equilibrium (1 = at equilibrium)
    Equilibrium_RPA = 0;    %test of whether RPA saturation is at equilibrium (1 = at equilibirium)
    while any([Equilibrium_RAD51,Equilibrium_RPA] == 0) == 1 & t(end) <= 25  %runs the whole time that the system is not at equilibrium
        Event_Count = Event_Count+1;    %advances event counter
        Event = Event_Count;
    %%% Search Process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Available_Counts,Bound_Counts,RAD51_Mon_SC,RAD51_Dim_SC,RAD51_Mon_I,RAD51_Dim_I,RAD51_Mon_DC,RAD51_Dim_DC,RPA_I,RPA_SC,RPA_DC,Free_for_RPA_D,Left_RAD51_Dimer_Filament] = LatticeSearch(N,n_RAD51,n_A,n_D,n_RPA,RPA_A,RPA_D,RAD51,RPA_D_HingedOpen,DNA);
        RAD51_Dim_BoundAtSpot(Left_RAD51_Dimer_Filament) = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Macro_Propensity = [k_on_RAD51*L_RAD51_M,k_on_RAD51*L_RAD51_M*w_RAD51,k_on_RAD51*L_RAD51_M*(w_RAD51^2),k_on_RAD51*L_RAD51_D,k_on_RAD51*L_RAD51_D*w_RAD51,k_on_RAD51*L_RAD51_D*(w_RAD51^2),k_on_RPA_A*L_RPA,k_on_RPA_A*L_RPA*w_RPA,k_on_RPA_A*(w_RPA^2)].*Available_Counts;  %propensity functions for macroscopic reactions
        Unbinding_Propensity = [k_off_RAD51,k_off_RAD51,k_off_RPA_A,k_off_RPA_D].*Bound_Counts;   %propensity functions of unbinding reactions
                    % 1 - RAD51 Monomer Unbinding
                    % 2 - RAD51 Dimer Unbinding
                    % 3 - RPA Macroscopically Dissociating
                    % 4 - RPA-D Microscopically Dissociating
        Micro_Propensity = k_on_RPA_D*numel(Free_for_RPA_D);   %propensity function for RPA-D rebinding (hinging closed)
        Full_Propensity = [Macro_Propensity, Unbinding_Propensity, Micro_Propensity];  %full listing of propensity functions for the simulation
                    % 1 - RAD51 Monomer I Binding
                    % 2 - RAD51 Monomer SC Binding
                    % 3 - RAD51 Monomer DC Binding
                    % 4 - RAD51 Dimer I Binding
                    % 5 - RAD51 Dimer SC Binding
                    % 6 - RAD51 Dimer DC Binding
                    % 7 - RPA Macro I Binding
                    % 8 - RPA Macro SC Binding
                    % 9 - RPA Macro DC Binding
                    % 10 - RAD51 Monomer Unbinding
                    % 11 - RAD51 Dimer Unbinding
                    % 12 - RPA Macro Dissociating
                    % 13 - RPA-D Micro Dissociating
                    % 14 - RPA-D Micro Binding

        a_0 = sum(Full_Propensity); %sum of all propensity functions
        Randoms = [rand,rand];  %random numbers for Monte Carlo step
        dt = (1/a_0)*log(1/Randoms(1)); %time until next reaction occurs

        %Direct Method of Calculating which reaction is next to occur
        if Full_Propensity(1) > Randoms(2)*a_0     %RAD51 Monomer I Binding Occurs
            j = 1;  %records which reaction occured
            RAD51_Mon_I_Bind_Spot = RAD51_Mon_I(randi(numel(RAD51_Mon_I)));    %selects random location to bind to that's available
            DNA(2,RAD51_Mon_I_Bind_Spot:RAD51_Mon_I_Bind_Spot+(n_RAD51-1)) = RAD51; %binds RAD51 monomer to DNA lattice
            RAD51_Mon_BoundAtSpot(RAD51_Mon_I_Bind_Spot) = 1;   %records that there is a protein bound at this location
            LocationHistory(1,Event_Count) = RAD51_Mon_I_Bind_Spot;   %stores location in LocationHistory
        elseif sum(Full_Propensity(1:2)) > Randoms(2)*a_0  %RAD51 Monomer SC Binding Occurs
            j = 2; %records which reaction occured
            RAD51_Mon_SC_Bind_Spot = RAD51_Mon_SC(randi(numel(RAD51_Mon_SC)));  %selects random location to bind RAD51 to that's available for corresponding reaction
            DNA(2,RAD51_Mon_SC_Bind_Spot:RAD51_Mon_SC_Bind_Spot+(n_RAD51-1)) = RAD51;   %binds RAD51 monomer
            RAD51_Mon_BoundAtSpot(RAD51_Mon_SC_Bind_Spot) = 1;  %records that a RAD51 Mon. is bound here
            LocationHistory(2,Event_Count) = RAD51_Mon_SC_Bind_Spot;    %records location of the event
        elseif sum(Full_Propensity(1:3)) > Randoms(2)*a_0  %RAD51 Monomer DC Binding Occurs
            j = 3; %records which reaction occured
            RAD51_Mon_DC_Bind_Spot = RAD51_Mon_DC(randi(numel(RAD51_Mon_DC)));  %selects random location for binding
            DNA(2,RAD51_Mon_DC_Bind_Spot:RAD51_Mon_DC_Bind_Spot+(n_RAD51-1)) = RAD51;   %binds RAD51 dimer to location
            RAD51_Mon_BoundAtSpot(RAD51_Mon_DC_Bind_Spot) = 1; %records where monomer of RAD51 is bound
            LocationHistory(3,Event_Count) = RAD51_Mon_DC_Bind_Spot;    %records where binding happened and when
        elseif sum(Full_Propensity(1:4)) > Randoms(2)*a_0  %RAD51 Dimer I Binding Occurs
            j = 4; %records which reaction occured
            RAD51_Dim_I_Bind_Spot = RAD51_Dim_I(randi(numel(RAD51_Dim_I))); %selects random location for dimer binding
            DNA(2,RAD51_Dim_I_Bind_Spot:RAD51_Dim_I_Bind_Spot+(2*n_RAD51-1)) = RAD51;   %binds RAD51 dimer
            RAD51_Mon_BoundAtSpot([RAD51_Dim_I_Bind_Spot,RAD51_Dim_I_Bind_Spot+n_RAD51]) = 1;   %records location of each monomer in the dimer
            LocationHistory(4,Event_Count) = RAD51_Dim_I_Bind_Spot; %records where location happened with this event
        elseif sum(Full_Propensity(1:5)) > Randoms(2)*a_0  %RAD51 Dimer SC Binding Occurs
            j = 5; %records which reaction occured
            RAD51_Dim_SC_Bind_Spot = RAD51_Dim_SC(randi(numel(RAD51_Dim_SC)));  %selects location for dimer binding
            DNA(2,RAD51_Dim_SC_Bind_Spot:RAD51_Dim_SC_Bind_Spot+(2*n_RAD51-1)) = RAD51; %binds dimer to location
            RAD51_Mon_BoundAtSpot([RAD51_Dim_SC_Bind_Spot,RAD51_Dim_SC_Bind_Spot+n_RAD51]) = 1; %reocrds location of each monomer part
            LocationHistory(5,Event_Count) = RAD51_Dim_SC_Bind_Spot;    %records where binding happened and which type of reaction
        elseif sum(Full_Propensity(1:6)) > Randoms(2)*a_0  %RAD51 Dimer DC Binding Occurs
            j = 6; %records which reaction occured
            RAD51_Dim_DC_Bind_Spot = RAD51_Dim_DC(randi(numel(RAD51_Dim_DC)));  %selects location for dimer binding
            DNA(2,RAD51_Dim_DC_Bind_Spot:RAD51_Dim_DC_Bind_Spot+(2*n_RAD51-1)) = RAD51; %binds dimer to location
            RAD51_Mon_BoundAtSpot([RAD51_Dim_DC_Bind_Spot,RAD51_Dim_DC_Bind_Spot+n_RAD51]) = 1; %reocrds location of each monomer part
            LocationHistory(6,Event_Count) = RAD51_Dim_DC_Bind_Spot;    %records where binding happened and which type of reaction
        elseif sum(Full_Propensity(1:7)) > Randoms(2)*a_0  %RPA Macro I Binding Occurs
            j = 7; %records which reaction occured
            RPA_Macro_I_Bind_Spot = RPA_I(randi(numel(RPA_I))); %chooses random location for RPA macro binding
            DNA(2,RPA_Macro_I_Bind_Spot:RPA_Macro_I_Bind_Spot+(n_A-1)) = RPA_A; %binds A part of RPA
            DNA(2,RPA_Macro_I_Bind_Spot+n_A:RPA_Macro_I_Bind_Spot+n_A+(n_D-1)) = RPA_D; %binds D component of RPA to lattice
            RPA_A_BoundAtSpot(RPA_Macro_I_Bind_Spot) = 1;   %records where A component is bound
            RPA_D_BoundAtSpot(RPA_Macro_I_Bind_Spot+n_A) = 1;  %records wehre D compoenet of RPA is bound
            LocationHistory(7,Event_Count) = RPA_Macro_I_Bind_Spot; %records where this reaction occured in this event
        elseif sum(Full_Propensity(1:8)) > Randoms(2)*a_0  %RPA Macro SC Binding Occurs
            j = 8; %records which reaction occured
            RPA_Macro_SC_Bind_Spot = RPA_SC(randi(numel(RPA_SC))); %chooses random location for RPA macro binding
            DNA(2,RPA_Macro_SC_Bind_Spot:RPA_Macro_SC_Bind_Spot+(n_A-1)) = RPA_A; %binds A part of RPA
            DNA(2,RPA_Macro_SC_Bind_Spot+n_A:RPA_Macro_SC_Bind_Spot+n_A+(n_D-1)) = RPA_D; %binds D component of RPA to lattice
            RPA_A_BoundAtSpot(RPA_Macro_SC_Bind_Spot) = 1;   %records where A component is bound
            RPA_D_BoundAtSpot(RPA_Macro_SC_Bind_Spot+n_A) = 1;  %records wehre D compoenet of RPA is bound
            LocationHistory(8,Event_Count) = RPA_Macro_SC_Bind_Spot; %records where this reaction occured in this event
        elseif sum(Full_Propensity(1:9)) > Randoms(2)*a_0  %RPA Macro DC Binding Occurs
            j = 9; %records which reaction occured
            RPA_Macro_DC_Bind_Spot = RPA_DC(randi(numel(RPA_DC))); %chooses random location for RPA macro binding
            DNA(2,RPA_Macro_DC_Bind_Spot:RPA_Macro_DC_Bind_Spot+(n_A-1)) = RPA_A; %binds A part of RPA
            DNA(2,RPA_Macro_DC_Bind_Spot+n_A:RPA_Macro_DC_Bind_Spot+n_A+(n_D-1)) = RPA_D; %binds D component of RPA to lattice
            RPA_A_BoundAtSpot(RPA_Macro_DC_Bind_Spot) = 1;   %records where A component is bound
            RPA_D_BoundAtSpot(RPA_Macro_DC_Bind_Spot+n_A) = 1;  %records wehre D compoenet of RPA is bound
            LocationHistory(9,Event_Count) = RPA_Macro_DC_Bind_Spot; %records where this reaction occured in this event
        elseif sum(Full_Propensity(1:10)) > Randoms(2)*a_0 %RAD51 Monomer Unbinding Occurs
            j = 10; %records which reaction occured
            RAD51_Bound_Monomers = find(RAD51_Mon_BoundAtSpot == 1);    %list of all locations a RAD51 monomer is bound
            RAD51_M_Unbind_Spot = RAD51_Bound_Monomers(randi(numel(RAD51_Bound_Monomers))); %selects random monomer that will unbind
            DNA(2,RAD51_M_Unbind_Spot:RAD51_M_Unbind_Spot+(n_RAD51-1)) = 0; %unbinds the monomer
            RAD51_Mon_BoundAtSpot(RAD51_M_Unbind_Spot) = 0; %removes location from BountAtSpot
            LocationHistory(10,Event_Count) = RAD51_M_Unbind_Spot;  %records where this event occured
        elseif sum(Full_Propensity(1:11)) > Randoms(2)*a_0 %RAD51 Dimer Unbinding Occurs
            j = 11; %records which reaction occured
            RAD51_Bound_Dimers = Left_RAD51_Dimer_Filament; %list of all locations where a dimer is bound
            RAD51_D_Unbind_Spot = RAD51_Bound_Dimers(randi(numel(RAD51_Bound_Dimers))); %randomly selects dimer to unbind from lattice
            DNA(2,RAD51_D_Unbind_Spot:RAD51_D_Unbind_Spot+(2*n_RAD51-1)) = 0;   %unbinds the dimer
            RAD51_Dim_BoundAtSpot(RAD51_D_Unbind_Spot) = 0; %removes location from BoundAtSpot (not necessary)
            RAD51_Mon_BoundAtSpot([RAD51_D_Unbind_Spot,RAD51_D_Unbind_Spot+n_RAD51]) = 0; %removes location from RAD51 Mon. BoundAtSpot record array
            LocationHistory(11,Event_Count) = RAD51_D_Unbind_Spot;  %records where this event occured
        elseif sum(Full_Propensity(1:12)) > Randoms(2)*a_0 %RPA Macro Dissociation
            j = 12; %records which reaction occured
            RPA_Bound_Macro = find(RPA_A_BoundAtSpot == 1); %list of locations where RPA is macro. bound
            RPA_Macro_Unbind_Spot = RPA_Bound_Macro(randi(numel(RPA_Bound_Macro))); %selection of random RPA molecule to unbind
            DNA(2,RPA_Macro_Unbind_Spot:RPA_Macro_Unbind_Spot+(n_A-1)) = 0; %unbinds the A part of RPA
            if DNA(1,RPA_Macro_Unbind_Spot+n_A) == RPA_D  %if the D part of RPA is hinged open...
                DNA(1,RPA_Macro_Unbind_Spot+n_A:RPA_Macro_Unbind_Spot+n_A+n_D-1) = 0;   %...then clear that spot
                RPA_D_HingedOpen(RPA_Macro_Unbind_Spot+n_A) = 0;    %records that RPA is no longer hinged open here (it's not bound here at all)
            else
                DNA(2,RPA_Macro_Unbind_Spot+n_A:RPA_Macro_Unbind_Spot+n_A+n_D-1) = 0;   %otherwise clear the bound D part of RPA
                RPA_D_BoundAtSpot(RPA_Macro_Unbind_Spot+n_A) = 0; %clear location of RPA-D from the BoundAtSpot array
            end
            RPA_A_BoundAtSpot(RPA_Macro_Unbind_Spot) = 0;   %clear location of unbinding from RPA_A_BoundAtSpot
            LocationHistory(12,Event_Count) = RPA_Macro_Unbind_Spot;    %record where this event occured at
        elseif sum(Full_Propensity(1:13)) > Randoms(2)*a_0 %RPA-D Micro Dissociation
            j = 13; %records which reaction occured
            RPA_D_Bound = find(RPA_D_BoundAtSpot == 1); %all locations where RPA-D is bound to the lattice
            RPA_D_Micro_Unbind_Spot = RPA_D_Bound(randi(numel(RPA_D_Bound))); %random RPA-D to unbind (hinge open)
            DNA(2,RPA_D_Micro_Unbind_Spot:RPA_D_Micro_Unbind_Spot+(n_D-1)) = 0; %unbinds RPA-D from corresponding location
            DNA(1,RPA_D_Micro_Unbind_Spot:RPA_D_Micro_Unbind_Spot+(n_D-1)) = RPA_D; %stores RPA-D in top row of DNA to show that it's hinged open
            RPA_D_BoundAtSpot(RPA_D_Micro_Unbind_Spot) = 0; %clears the recording of that RPA-D part being bound there
            RPA_D_HingedOpen(RPA_D_Micro_Unbind_Spot) = 1;  %records that the RPA-D is hinged open here
            LocationHistory(13,Event_Count) = RPA_D_Micro_Unbind_Spot;  %records where this reaction occured
        elseif sum(Full_Propensity(1:14)) > Randoms(2)*a_0 %RPA-D Micro Re-binding
            j = 14; %records which reaction occured
            RPA_D_Micro_Bind_Spot = Free_for_RPA_D(randi(numel(Free_for_RPA_D)));  %random location where RPA-D can rebind (and will in this event)
            DNA(2,RPA_D_Micro_Bind_Spot:RPA_D_Micro_Bind_Spot+(n_D-1)) = RPA_D; %binds RPA-D to the DNA lattice
            DNA(1,RPA_D_Micro_Bind_Spot:RPA_D_Micro_Bind_Spot+(n_D-1)) = 0; %clears RPA-D from being hinged open
            RPA_D_BoundAtSpot(RPA_D_Micro_Bind_Spot) = 1;   %records that RPA-D is now bound to the DNA lattice
            RPA_D_HingedOpen(RPA_D_Micro_Bind_Spot) = 0;    %records that RPA-D is no longer hinged open
            LocationHistory(14,Event_Count) = RPA_D_Micro_Bind_Spot;    %records where this event occured
        end

        %Bug Checking - Finding Whole Proteins
        if numel(find(DNA(2,:) == RAD51))/n_RAD51 ~= round(numel(find(DNA(2,:) == RAD51))/n_RAD51)  %checks if there are an integer number of RAD51 proteins
            disp('ERROR: BROKEN RAD51');
            disp(['Last Reaction: ', num2str(j(end))]);
            break;
        elseif numel(find(DNA(2,:) == RPA_A))/n_A ~= round(numel(find(DNA(2,:) == RPA_A))/n_A)     %checks if there is an integer number of RPA-A bound
            disp('ERROR: BROKEN RPA-A');
            disp(['Last Reaction: ', num2str(j(end))]);
            break;
        elseif numel(find(DNA(2,:) == RPA_D))/n_D ~= round(numel(find(DNA(2,:) == RPA_D))/n_D)      %checks for an integer number of bound RPA_D
            disp('ERROR: BROKEN RPA-D (Bound)');
            disp(['Last Reaction: ', num2str(j(end))]);
            break;
        elseif numel(find(DNA(1,:) == RPA_D))/n_D ~= round(numel(find(DNA(1,:) == RPA_D))/n_D)      %checks for an integer number of hinged open RPA-D
            disp('ERROR: BROKEN RPA-D (Open)');
            disp(['Last Reaction: ', num2str(j(end))]);
            break;
        end

        t = [t,t(end)+dt];  %advance time according to time interval selected by Gillespie Algorithm

        FracCover_RAD51 = [FracCover_RAD51,numel(find(DNA(2,:) == RAD51))/N]; %RAD51 saturation of the DNA lattice after each event
        FracCover_RPA_A = [FracCover_RPA_A,numel(find(DNA(2,:) == RPA_A))/N];  %saturation of RPA-A
        FracCover_RPA_D = [FracCover_RPA_D,numel(find(DNA(2,:) == RPA_D))/N];  %saturation of RPA-D on the lattice
        FracCover_RPA = [FracCover_RPA,sum([FracCover_RPA_A(Event_Count+1),FracCover_RPA_D(Event_Count+1)])];  %saturation of all parts of RPA
        FracCover_Total = [FracCover_Total,sum([FracCover_RAD51(Event_Count+1),FracCover_RPA(Event_Count+1)])];  %total saturation of the lattice (both RPA and RAD51)

    % Equilibrium Testing - Linear Slope & Intercept Method %%%%%%%%%%%%%%%%%%%
        if Event_Count >= 1000 %only tests for equilibrium after 1000 events have occured
            t_Equilibrium_Test = t(end-round(0.25*(Event_Count+1)):end);  %time values that we're testing for equilibrium
            RPA_Equilibrium_Test = FracCover_RPA((end-round(0.25*(Event_Count+1))):end);   %last 1/4 of Events saturation data for RPA
            RAD51_Equilibrium_Test = FracCover_RAD51((end-round(0.25*(Event_Count+1))):end);   %last 1/4 of Events saturation data for RAD51

            RPA_Avg_Saturation_Holder = sum(RPA_Equilibrium_Test)/numel(RPA_Equilibrium_Test); %average saturation in last 1/4 of Events (RPA)
            RAD51_Avg_Saturation_Holder = sum(RAD51_Equilibrium_Test)/numel(RAD51_Equilibrium_Test);   %average saturation in last 1/4 of Events (RAD51)

            RAD51_Fit = polyfit(t_Equilibrium_Test,RAD51_Equilibrium_Test,1);   %linear fit for RAD51 data (slope, y-int)
            RAD51_Yint_Error = abs(RAD51_Avg_Saturation_Holder-RAD51_Fit(2))/RAD51_Avg_Saturation_Holder; %y-intercept error of linear fit for RAD51 data
            RPA_Fit = polyfit(t_Equilibrium_Test,RPA_Equilibrium_Test,1);   %linear fit to RPA data (slope, y-int)
            RPA_Yint_Error = abs(RPA_Avg_Saturation_Holder-RPA_Fit(2))/RPA_Avg_Saturation_Holder;    %y-intercept error compared to average RPA saturation

            if abs(RPA_Fit(1)) < 0.01 & (RPA_Yint_Error < 0.05 | isnan(RPA_Yint_Error))   %if slope of RPA data is essentially zero and y-intercept is very close to avg. saturation value... (slope limit is change in saturation of 1% (~17 proteins) per 1 time interval)
                Equilibrium_RPA = 1;    %...then at equilibrium
            else
                Equilibrium_RPA = 0;    %...otherwise reset to not at equilibrium
            end
            if abs(RAD51_Fit(1)) < 0.01 & (RAD51_Yint_Error < 0.05 | isnan(RAD51_Yint_Error)) %if the slope of RAD51 data is essentially zero and y-intercept is very close to avg. saturation value... (slope limit is change in saturation of 1% (~3 proteins) per 1 time interval)
                Equilibrium_RAD51 = 1;  %...then we're at equilibrium
            else
                Equilibrium_RAD51 = 0;    %...otherwise reset to not at equilibrium
            end
        end   
    end         %end of a single simulation
    
    RPA_Avg_Saturation(Simulations) = RPA_Avg_Saturation_Holder;    %records equilibrium avg. RPA saturation
    RAD51_Avg_Saturation(Simulations) = RAD51_Avg_Saturation_Holder;    %records equilibrium avg. RAD51 saturation
    Total_Avg_Saturation(Simulations) = sum(FracCover_Total(end-round(0.25*Event_Count):end))/numel(FracCover_Total(end-round(0.25*Event_Count):end));  %records equilibrium avg. for the whole DNA molecule
    t_Equilibrium(Simulations) = t(end-round(0.25*(Event_Count+1)));   %time where equilibrium occured

%     Record data for growth profiles
    time_Data{Simulations} = t;   %records time data for each simulation
    FracCover_RAD51_Data{Simulations} = FracCover_RAD51;    %records RAD51 saturation throughout each simulation
    FracCover_RPA_A_Data{Simulations} = FracCover_RPA_A;    %records RPA-A saturation throughout each simulation
    FracCover_RPA_D_Data{Simulations} = FracCover_RPA_D;    %reocrds RPA-D saturation throughout each simulation
    FracCover_RPA_Data{Simulations} = FracCover_RPA;    %records total RPA saturation throughout each simulation
    FracCover_Total_Data{Simulations} = FracCover_Total;    %records the total saturation throughout each simulation

    Simulation_Times(Simulations) = toc;    %times each simulation
    send(q_Count,Simulations); %sends simulation count to q_Count and updates progress bar
end
delete(w);  %closes waitbar

[X_Heatmap,Y_Heatmap] = meshgrid(L_RAD51_Total_Values,k_off_RPA_A_Values); %generates X and Y values for heatmaps
Z_Heatmap_RAD51 = griddata(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),RAD51_Avg_Saturation,X_Heatmap(:),Y_Heatmap(:));  %generates heatmap data based on RAD51 saturation
Z_Heatmap_RAD51 = reshape(Z_Heatmap_RAD51,size(X_Heatmap)); %reshapes RAD51 saturation data into appropriately sized matrix
Z_Heatmap_RPA = griddata(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),RPA_Avg_Saturation,X_Heatmap(:),Y_Heatmap(:));  %generates heatmap data based on RPA saturation
Z_Heatmap_RPA = reshape(Z_Heatmap_RPA,size(X_Heatmap)); %reshapes RPA saturation data into appropriately sized matrix
Z_Heatmap_Total = griddata(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),Total_Avg_Saturation,X_Heatmap(:),Y_Heatmap(:));  %generates heatmap data based on the total saturation (sum of RPA and RAD51)
Z_Heatmap_Total = reshape(Z_Heatmap_Total,size(X_Heatmap)); %reshapes Total saturation data into appropriately sized matrix
Z_Heatmap_t = griddata(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),t_Equilibrium,X_Heatmap(:),Y_Heatmap(:));    %generates heatmap data based on the time to equilibrium
Z_Heatmap_t = reshape(Z_Heatmap_t,size(X_Heatmap)); %reshapes time to equilibrium data into the appropriately sized matrix

figure;
subplot(2,2,1); %RAD51 Saturation Heatmap
imagesc(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),Z_Heatmap_RAD51);
title('RAD51 Saturation');
xlabel('L_R_A_D_5_1 (Total)');
ylabel('RPA-A k_o_f_f');
colorbar;
box on;
subplot(2,2,2); %RPA Saturation Heatmap
imagesc(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),Z_Heatmap_RPA);
title('RPA Saturation');
xlabel('L_R_A_D_5_1 (Total)');
ylabel('RPA-A k_o_f_f');
colorbar;
box on;
subplot(2,2,3); %Total Saturation Heatmap
imagesc(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),Z_Heatmap_Total);
title('Total Saturation');
xlabel('L_R_A_D_5_1 (Total)');
ylabel('RPA-A k_o_f_f');
colorbar;
box on;
subplot(2,2,4); %Time to Equilibrium Heatmap
imagesc(Heatmap_Parameters(1,:),Heatmap_Parameters(2,:),Z_Heatmap_t);
title('Time to Equilibrium');
xlabel('L_R_A_D_5_1 (Total)');
ylabel('RPA-A k_o_f_f');
colorbar;
box on;

toc