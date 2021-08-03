function [Available_Counts,Bound_Counts,RAD51_Mon_SC,RAD51_Dim_SC,RAD51_Mon_I,RAD51_Dim_I,RAD51_Mon_DC,RAD51_Dim_DC,RPA_I,RPA_SC,RPA_DC,Free_for_RPA_D,Left_RAD51_Dimer_Filament] = LatticeSearch(N,n_RAD51,n_A,n_D,n_RPA,RPA_A,RPA_D,RAD51,RPA_D_HingedOpen,DNA)
%Function to search the DNA lattice for all populations and counts of
%available and bount proteins. Runs with each Event in a simulation. Main
%use is so that it can be adjusted and applied to all codes (via GitHub)
%simultaneously and easily. Copy/Paste of what the code was in the script.

    Gap_Left = find(diff([1 DNA(2,:) 1])<0 & diff([1 DNA(2,:) 1]) ~= RPA_A-RPA_D & diff([1 DNA(2,:) 1]) ~=RPA_D-RAD51 & diff([1 DNA(2,:) 1]) ~= RPA_A-RAD51);    %left most available location of all gaps on lattice
    Gap_Right = find(diff([1 DNA(2,:) 1])>0 & diff([1 DNA(2,:) 1]) ~= RPA_D-RPA_A & diff([1 DNA(2,:) 1]) ~= RAD51-RPA_D & diff([1 DNA(2,:) 1]) ~= RAD51-RPA_A)-1; %right most available location of all gaps on lattice
    Gap_Size = (Gap_Right-Gap_Left)+1;    %calculates the size of each gap
    Gap_Edges = [Gap_Left;Gap_Right];   %lists all gap edges (1st row = left edge, 2nd row = right edge)

    RAD51_M_Available_Gap_Edges = [Gap_Left(Gap_Size >= n_RAD51);Gap_Right(Gap_Size >= n_RAD51)];   %lists the left (1st row) and right (2nd row) edges
                                                                                                    %of gaps large enough for RAD51 monomers
    RAD51_D_Available_Gap_Edges = [Gap_Left(Gap_Size >= 2*n_RAD51);Gap_Right(Gap_Size >= 2*n_RAD51)];   %lists the left (1st row) and right (2nd row) edges
                                                                                                        %of gaps large enough for RAD51 dimers
    RPA_Available_Gap_Edges = [Gap_Left(Gap_Size >= n_RPA);Gap_Right(Gap_Size >= n_RPA)];   %lists left (1st row) and right (2nd row) edges of gaps
                                                                                            %large enough for an RPA protein to bind

    %Doubly Contiguous Search - Easiest
    RAD51_Mon_DC = Gap_Left(Gap_Size == n_RAD51 & Gap_Left > 1 & Gap_Left < N-(n_RAD51-1));   %all available RAD51 Mon. DC sites
    RAD51_Dim_DC = Gap_Left(Gap_Size == 2*n_RAD51 & Gap_Left > 1 & Gap_Left < N-(2*n_RAD51-1)); %all available RAD51 Dim. DC sites
    RPA_DC = Gap_Left(Gap_Size == n_RPA & Gap_Left > 1 & Gap_Left < N-(n_RPA-1));   %all available RPA DC sites

    %Singly Contiguous Search
    RAD51_Mon_SC = unique([RAD51_M_Available_Gap_Edges(1,:),RAD51_M_Available_Gap_Edges(2,:)-(n_RAD51-1)]); %position of all locations at the edges of gaps available for RAD51 Mon.
    RAD51_Mon_SC(ismember(RAD51_Mon_SC,RAD51_Mon_DC)) = [];     %clears location if it's already been recorded as a DC location
    RAD51_Mon_SC(ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == RPA_A)-n_RAD51) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == RPA_D)-n_RAD51) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == -RPA_A)) == 1 | ismember(RAD51_Mon_SC,find(diff([0 DNA(2,:) 0]) == -RPA_D)) == 1) = []; %clears any position that neighbors a protein other than RAD51 
    if DNA(2,n_RAD51+1) ~= RAD51    %if the first location isn't SC...
        RAD51_Mon_SC(RAD51_Mon_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-n_RAD51) ~= RAD51       %if the final available location on lattice isn't SC...
        RAD51_Mon_SC(RAD51_Mon_SC == N-(n_RAD51-1)) = [];   %...clear it
    end
    RAD51_Dim_SC = unique([RAD51_D_Available_Gap_Edges(1,:),RAD51_D_Available_Gap_Edges(2,:)-(2*n_RAD51-1)]);   %positions of all locations at the edges of gaps available for RAD51 Dim.
    RAD51_Dim_SC(ismember(RAD51_Dim_SC,RAD51_Dim_DC)) = [];     %clears locations that have already been counted as a DC location
    RAD51_Dim_SC(ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == RPA_A)-(2*n_RAD51)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == RPA_D)-(2*n_RAD51)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == -RPA_A)) == 1 | ismember(RAD51_Dim_SC,find(diff([0 DNA(2,:) 0]) == -RPA_D)) == 1) = []; %clears all position that neighbor a protein other than RAD51
    if DNA(2,2*n_RAD51+1) ~= RAD51    %if the first location isn't SC...
        RAD51_Dim_SC(RAD51_Dim_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-(2*n_RAD51)) ~= RAD51       %if the final available location on lattice isn't SC...
        RAD51_Dim_SC(RAD51_Dim_SC == N-(2*n_RAD51-1)) = [];   %...clear it
    end
    RPA_SC = unique([RPA_Available_Gap_Edges(1,:),RPA_Available_Gap_Edges(2,:)-(n_RPA-1)]);   %positions of all locations at the edges of gaps available for RPA
    RPA_SC(ismember(RPA_SC,RPA_DC)) = [];     %clears locations that have already been counted as a DC location
    RPA_SC(ismember(RPA_SC,find(diff([0 DNA(2,:) 0]) == RAD51)-n_RPA) == 1 | ismember(RPA_SC,find(diff([0 DNA(2,:) 0]) == -RAD51)) == 1) = []; %clears all locations that don't neighbor an RPA protein (A or D)
    if DNA(2,n_RPA+1) ~= RPA_A   %if the first location isn't SC...
        RPA_SC(RPA_SC == 1) = [];   %...clear it
    end
    if DNA(2,N-n_RPA) ~= RPA_D & DNA(2,N-n_RPA) ~= RPA_A     %if the final available location on lattice isn't SC...
        RPA_SC(RPA_SC == N-(n_RPA-1)) = [];   %...clear it
    end

    %Isolated Search
    RAD51_Mon_I = [];   %initializes an empty array for RAD51 monomer isolated available sites
    Gap_Size_RAD51_M_I = Gap_Size(Gap_Size > n_RAD51+1);    %gap sizes of gaps large enough for a RAD51 monomer to bind in an isolated position
    RAD51_M_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RAD51_M_I) == 1);  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for a = 1:length(RAD51_M_Available_I_Gap_Edges(1,:))
        RAD51_Mon_I = [RAD51_Mon_I, RAD51_M_Available_I_Gap_Edges(1,a)+1:1:RAD51_M_Available_I_Gap_Edges(2,a)-n_RAD51]; %adds sequential locations that are isolated locations available for RAD51 Mon
    end
    if DNA(2,n_RAD51) ~= RAD51 & DNA(2,1:n_RAD51) == 0   %if left end of lattice is available and isn't SC...
        RAD51_Mon_I = [1,RAD51_Mon_I];  %...add the first location to the list
    end
    if DNA(2,N-n_RAD51) ~= RAD51 & DNA(2,N-(n_RAD51-1):end) == 0   %if the right end of the lattice is available and isn't SC...
        RAD51_Mon_I = [RAD51_Mon_I,N-(n_RAD51-1)];  %...add the last location to the list
    end
    RAD51_Dim_I = [];   %initializes an empty array for RAD51 dimer isolated available sites
    Gap_Size_RAD51_D_I = Gap_Size(Gap_Size > 2*n_RAD51+1);    %gap sizes large enough for RAD51 Dimers
    RAD51_D_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RAD51_D_I));  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for b = 1:length(RAD51_D_Available_I_Gap_Edges(1,:))
        RAD51_Dim_I = [RAD51_Dim_I, RAD51_D_Available_I_Gap_Edges(1,b)+1:1:RAD51_D_Available_I_Gap_Edges(2,b)-(2*n_RAD51)]; %adds sequential locations that are isolated locations available for RAD51 Dim
    end
    if DNA(2,1:2*n_RAD51) == 0 & DNA(2,(2*n_RAD51)+1) ~= RAD51   %if left end of lattice is available and isn't SC...
        RAD51_Dim_I = [1,RAD51_Dim_I];  %...add the first location to the list
    end
    if DNA(2,N-(2*n_RAD51-1):end) == 0 & DNA(2,N-(2*n_RAD51)) ~= RAD51   %if the right end of the lattice isn't SC...
        RAD51_Dim_I = [RAD51_Dim_I,N-(2*n_RAD51-1)];  %...add the last location to the list
    end
    RPA_I = [];   %initializes an empty array for RPA isolated available sites
    Gap_Size_RPA = Gap_Size(Gap_Size > n_RPA);  %gap sizes which are large enough for RPA isolated binding
    RPA_Available_I_Gap_Edges = Gap_Edges(:,ismember(Gap_Size,Gap_Size_RPA));  %lists of gap edges (left = 1st row, right = 2nd row) that are large enough for any isolated binding
    for c = 1:length(RPA_Available_I_Gap_Edges(1,:))
        RPA_I = [RPA_I, RPA_Available_I_Gap_Edges(1,c)+1:1:RPA_Available_I_Gap_Edges(2,c)-n_RPA]; %adds sequential locations that are isolated locations available for RAD51 Mon
    end
    if DNA(2,1:n_RPA) == 0 & DNA(2,n_RPA+1) ~= RPA_A   %if left end of lattice is available and isn't SC...
        RPA_I = [1,RPA_I];  %...add the first location to the list
    end
    if DNA(2,N-(n_RPA-1):end) == 0 & DNA(2,N-(n_RPA+1)) ~= RPA_A & DNA(2,N-(n_RPA+1)) ~= RPA_D  %if the right end of the lattice is available and isn't SC...
        RPA_I = [RPA_I,N-(n_RPA-1)];  %...add the last location to the list
    end

    %RPA-D Availble Locations
    Free_for_RPA_D = [];    %initializes an array for where hinged open RPA-D can rebind to DNA
    for d = find(RPA_D_HingedOpen == 1)   %check each location where RPA-D is hinged open
        if DNA(2,d:d+(n_D-1)) == 0     %if no proteins are bound below hinged open RPA-D...
            Free_for_RPA_D = [Free_for_RPA_D,d];    %...adds location to store where RPA-D can rebind this event
        end
    end
    %Population of Bound Proteins
    x_Bound_RAD51_M = numel(find(DNA(2,:) == RAD51))/n_RAD51;    %calculates how many RAD51 Monomers are actively bound to DNA
    x_Bound_RPA_A = numel(find(DNA(2,:) == RPA_A))/n_A;    %calculates the number of bound RPA-A
    x_Bound_RPA_D = numel(find(DNA(2,:) == RPA_D))/n_D;    %calculates the number of actively bound RPA-D
    %Search for RAD51 Dimers (Uses Filament Lengths)
    RAD51_Dim_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Dimers are bound on the lattice
    RAD51_Filament_Edges = [find(diff([0 DNA(2,:) 0]) == 51 | diff([0 DNA(2,:) 0]) == 50 | diff([0 DNA(2,:) 0]) == 48);find(diff([0 DNA(2,:) 0]) == -51 | diff([0 DNA(2,:) 0]) == -50 | diff([0 DNA(2,:) 0]) == -48)-1];  %edges of RAD51 filaments (left = 1st row; right = 2nd row)
    RAD51_Filament_Lengths = RAD51_Filament_Edges(2,:)-RAD51_Filament_Edges(1,:)+1; %length of each RAD51 Filament on the lattice
    RAD51_D_Filament_Locations = RAD51_Filament_Edges(1,RAD51_Filament_Lengths > n_RAD51);  %location of all filaments which contain a dimer
    RAD51_D_Filament_Lengths = RAD51_Filament_Lengths(RAD51_Filament_Lengths > n_RAD51);   %length of RAD51 filaments which contain a dimer
    Monomers_per_Dimer_Filament = RAD51_D_Filament_Lengths./n_RAD51;    %number of monomers in each filament that contains a dimer
    Left_RAD51_Dimer_Filament = []; %initializes array to record locations of dimers in filaments
    if numel(RAD51_D_Filament_Locations) > 0   %if there are any filaments containing dimers...
        for e = 1:numel(RAD51_D_Filament_Locations)     %...check each one to add it to the list of dimer locations
            if RAD51_D_Filament_Lengths(e) == 2*n_RAD51    %if there's only one dimer in the filament...
                Left_RAD51_Dimer_Filament = [Left_RAD51_Dimer_Filament,RAD51_D_Filament_Locations(e)];  %...then record the location of that dimer   %...then record that there is a dimer bound at corresponding location
            else   %...otherwise we have to check each monomer location
                for f = 1:Monomers_per_Dimer_Filament(e)-1     %record that there's a possible dimer at each monomer location (except the last one)
                    Left_RAD51_Dimer_Filament = [Left_RAD51_Dimer_Filament,RAD51_D_Filament_Locations(e)+((f-1)*n_RAD51)];  %record location of the beginning of each dimer within the filament
                end
            end
        end
    end
    RAD51_Dim_BoundAtSpot(Left_RAD51_Dimer_Filament) = 1;   %records where all possible dimers are located
    x_Bound_RAD51_D = numel(find(RAD51_Dim_BoundAtSpot == 1)); %number of RAD51 Dimers bound to lattice
    
    Available_Counts = [numel(RAD51_Mon_I),numel(RAD51_Mon_SC),numel(RAD51_Mon_DC),numel(RAD51_Dim_I),numel(RAD51_Dim_SC),numel(RAD51_Dim_DC),numel(RPA_I),numel(RPA_SC),numel(RPA_DC)];
                %populations of the available locations listed in order
                % 1 - RAD51 Monomer Isolated
                % 2 - RAD51 Monomer Singly Contiguous
                % 3 - RAD51 Monomer Doubly Contiguous
                % 4 - RAD51 Dimer Isolated
                % 5 - RAD51 Dimer Singly Contiguous
                % 6 - RAD51 Dimer Doubly Contiguous
                % 7 - RPA Isolated
                % 8 - RPA Singly Contiguous
                % 9 - RPA Doubly Contiguous
     Bound_Counts = [x_Bound_RAD51_M,x_Bound_RAD51_D,x_Bound_RPA_A,x_Bound_RPA_D];  %counts of how many bound proteins of each type
                % 1 - RAD51 Monomer
                % 2 - RAD51 Dimer
                % 3 - RPA-A
                % 4 - RPA-D
                
end

