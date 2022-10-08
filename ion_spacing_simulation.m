% 0. CLEARING THE COMMAND WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

% 1. SETTING SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1.1 Setting Control Ion (Phosphorus)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % State and Radius (Units: Angstroms)
        controlIonIsExcited = true;
        if (controlIonIsExcited)
            controlRadius = 100;            % Estimation of the control qubit's radius in its excited state
        else
            controlRadius = 1.64;           % The wavefunction radius of the control qubit in its excited state
        end
    
        % Energy (Units: keV)
        controlIonEnergy = 6;
        
    % 1.2 Identifying Target Depth to Set Non Control Ion Energies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p6 = (4.0745059722)*(10^(-14)); p5 = (-3.3601039663)*(10^(-11)); p4 = (-7.3049769110)*(10^(-9)); p3 = (1.1186615874)*(10^(-5)); p2 = (-4.2105628811)*(10^(-3)); p1 = (1.3398417218)*(10^(1)); p0 = (3.2060314440)*(10^(1));
        depth_function_P = @(x) p6 * x.^6 + p5 * x.^5 + p4 * x.^4 + p3 * x.^3 + p2 * x.^2 + p1 * x + p0;
        target_depth = depth_function_P(controlIonEnergy);

    % 1.3 Setting Non-Control Ion 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Type ('P' or 'Sb' or 'Bi')
        nonControlIonType = 'Bi';
        
        % Energy & Radius (Units: keV , Angstroms)
        if (nonControlIonType == 'P')
            nonControlRadius = 1.64;
            nonControlIonEnergy = controlIonEnergy;
        elseif (nonControlIonType == 'Sb')
            nonControlRadius = 2.33;
            sb6 = (-8.8565105611)*(10^(-13)); sb5 = (1.4252712448)*(10^(-9)); sb4 = (-8.8606740568)*(10^(-7)); sb3 = (2.6883211211)*(10^(-4)); sb2 = (-4.1682800162)*(10^(-2)); sb1 = (7.1589482290)*(10^(0)); sb0 = (4.7729568585)*(10^(1));
            depth_function_Sb = @(x) sb6 * x.^6 + sb5 * x.^5 + sb4 * x.^4 + sb3 * x.^3 + sb2 * x.^2 + sb1 * x + sb0 - target_depth;
            nonControlIonEnergy = fzero(depth_function_Sb , 0);
        elseif (nonControlIonType == 'Bi')
            nonControlRadius = 2.8;
            bi6 = (-9.6413455226)*(10^(-13)); bi5 = (1.5479457319)*(10^(-9)); bi4 = (-9.6481992870)*(10^(-7)); bi3 = (2.9614565783)*(10^(-4)); bi2 = (-4.7774657090)*(10^(-2)); bi1 = (6.4514781640)*(10^(0)); bi0 = (5.6604602995)*(10^(1));   
            depth_function_Bi = @(x) bi6 * x.^6 + bi5 * x.^5 + bi4 * x.^4 + bi3 * x.^3 + bi2 * x.^2 + bi1 * x + bi0 - target_depth;
            nonControlIonEnergy = fzero(depth_function_Bi , 0);
        end


    % 1.4 Setting Final Parameters 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The range of spacings tested in the simulation and the number of arrays tested per spacing value
        minSpacing = 20; maxSpacing = 250; spacingStep = 1;
        spacingsRange = [minSpacing:spacingStep:maxSpacing];
        arrays_per_spacing = 1000;
    
        % Overlaps per spacing
        overlap_Average_per_Spacing = zeros(length(spacingsRange),1);
        overlap_Disparity_per_Spacing = zeros(length(spacingsRange),1);
        overlap_Ideal_per_Spacing = zeros(length(spacingsRange),1);

% 2. SIMULATING ALL THE DIFFERENT SPACINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for currentSpacing = spacingsRange

    % 2.1 Initialising current spacing parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spacing ID
        spacingID = int16(((currentSpacing - minSpacing) + spacingStep) * (1/spacingStep));
        test_text = 'Current Spacing ID = ' + string(spacingID);
        disp(test_text);

        % Impact Locations for current spacing
        ImpactLocations = [0 1 1 ; 0 2 1 ; 0 3 1 ; 0 1 2; 0 2 2; 0 3 2; 0 1 3; 0 2 3; 0 3 3]; 
        ImpactLocations = ImpactLocations.*currentSpacing;

        % Overlap per array in current spacing
        overlap_per_Array = zeros(arrays_per_spacing,1);

        % Ideal Overlap
        overlap_Ideal = 0;
        for ion1 = 1 : 9
            for ion2 = 1 : 9
                if (ion1 <  ion2)
                    overlap_Ideal_per_Spacing(spacingID) = overlap_Ideal_per_Spacing(spacingID) + ( spherical_overlap(ImpactLocations , ion1 , ion2, nonControlRadius, controlRadius , 0) / 36 );
                end
            end
        end

    % 2.2 Simulating Arrays AND Calculating Overlaps (Average & Disparity) in Current Spacing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for currentArray = 1 : arrays_per_spacing
        
            % 2.2.1 Initialising current array parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            overlap_IonPairs_per_Array = zeros(36,1);
            DestinationLocations = ImpactLocations;
    
            % 2.2.2 Implantating Individual 3x3 Arrays
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Initialising control condition
                is_control_qubit = false;
                
                % Iterating through ions
                for currentIon = 1 : 9
                    % Setting middle ion as control ion
                    if (currentIon == 5)
                        is_control_qubit = true;
                    end
                        
                    % Calculating current ion paths
                    Xd = path_along_axis('x' , nonControlIonType , nonControlIonEnergy , controlIonEnergy , is_control_qubit);
                    Yd = path_along_axis('y' , nonControlIonType , nonControlIonEnergy , controlIonEnergy , is_control_qubit);
                    Zd = path_along_axis('z' , nonControlIonType , nonControlIonEnergy , controlIonEnergy , is_control_qubit);
                    
                    % Starting Location + Path = Destination
                    DestinationLocations(currentIon,1) = ImpactLocations(currentIon,1) + Xd;
                    DestinationLocations(currentIon,2) = ImpactLocations(currentIon,2) + Yd;
                    DestinationLocations(currentIon,3) = ImpactLocations(currentIon,3) + Zd;
                    
                    % Resetting control condition
                    is_control_qubit = false;
                end
        
            % 2.2.3 Calculating Overlaps per Array
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ionPairID = 1;
                for ion1 = 1 : 9
                   for ion2 = 1 : 9
                       if (ion1 < ion2)
                           overlap_IonPairs_per_Array(ionPairID) = spherical_overlap(DestinationLocations , ion1 , ion2 , nonControlRadius , controlRadius , ionPairID);
                           overlap_per_Array(currentArray) = overlap_per_Array(currentArray) + ( overlap_IonPairs_per_Array(ionPairID) / 36);
                           ionPairID = ionPairID + 1;
                       end
                   end
                end
                
            % 2.2.4 Calculating Overlaps per Spacing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                overlap_Average_per_Spacing(spacingID) = overlap_Average_per_Spacing(spacingID) + ( overlap_per_Array(currentArray) / arrays_per_spacing);
                overlap_Disparity_per_Spacing(spacingID) = overlap_Ideal_per_Spacing(spacingID) - overlap_Average_per_Spacing(spacingID);
    
            % 2.2.5 Printing current array trajectories 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %print_tragectories(nonControlIonType , nonControlIonEnergy , controlIonEnergy , currentSpacing , controlRadius , nonControlRadius , ImpactLocations , DestinationLocations);
    
            % 2.2.6 Probing Progress
            %%%%%%%%%%%%%%%%%%%%%%%%
                if ( mod(currentArray,10) == 0 )
                    progressProbe = 'SpacingID: ' + string(spacingID) + '/' + string (length(spacingsRange)) +'  Spacing Value = ' + string(currentSpacing) + 'Å  Completion: ' + string( currentArray*100 / arrays_per_spacing ) + '%';
                    disp(progressProbe);
                end
        end    
end


% 3. OVERLAP PLOTS
%%%%%%%%%%%%%%%%%%
figure
hold on

    % 3.3 Initialising Figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_title = "Ideal Overlap";
    title(this_title,'FontWeight','bold','FontSize',18);
    xAxisTitle = "Spacing (Å)";
    yAxisTitle = "Volume (Å^3)";

    % 3.4 Plotting figure
    %%%%%%%%%%%%%%%%%%%%%%%
        xlabel(xAxisTitle, 'color' , 'black' ,'FontSize',14);
        ylabel(yAxisTitle, 'color' , 'black' ,'FontSize',14);

        scatter(spacingsRange , overlap_Average_per_Spacing , 'blue' , 'o' , 'filled' );
        scatter(spacingsRange , overlap_Disparity_per_Spacing , 'red' , 'o' , 'filled' );
        scatter(spacingsRange , overlap_Ideal_per_Spacing , 'green' , 'o' , 'filled' );
    
        lgd = legend('Average Overlap','Overlap Disparity','Ideal Overlap','location','southoutside');
        lgd.FontSize = 12;
        lgd.NumColumns = 3;

hold off
