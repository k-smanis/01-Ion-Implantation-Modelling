%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The nonControlIonEnergy variable is only useful in that the function's prototype requires an input for it. Other than that it has no use.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. CLEARING THE COMMAND WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

% 1. INITIALISING SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
controlRadius = 100; is_control_qubit = true; nonControlIonEnergy = 0;
minEnergy = 2; maxEnergy = 30; energyStep = 0.2; energyLevels = [minEnergy:energyStep:maxEnergy];
acceptabilityRates = zeros(1,length(energyLevels));
averageImplantationDepths = zeros(1,length(energyLevels));
iterations_per_loop = 500000;

% 2. RUNNING VERYING ENERGY LOOPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for controlIonEnergy = energyLevels

    % 2.1 Identifying the energy loop (ranges from 1 to 17)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energyLoopID = int16(((controlIonEnergy - minEnergy) + energyStep) * (1/energyStep));

    % 2.2 Calculating the Average Depth of the specific energy level
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r1 = (4.0745059722)*(10^(-14)); r2 = (-3.3601039663)*(10^(-11)); r3 = (-7.3049769110)*(10^(-9)); r4 = (1.1186615874)*(10^(-5)); r5 = (-4.2105628811)*(10^(-3)); r6 = (1.3398417218)*(10^(1)); range_const = (3.2060314440)*(10^(1));
    range_function = @(x) r1 * x^6 + r2 * x^5 + r3 * x^4 + r4 * x^3 + r5 * x^2 + r6 * x + range_const;
    averageImplantationDepths(energyLoopID) = range_function(controlIonEnergy);  

    % 2.3 Energy loop of 1000 iterations starts here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for currentIteration = 1:iterations_per_loop
    
        % 2.3.1 Simulating single P ion implantation (in terms of depth only)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        depth = path_along_axis('x' , 'Phosphorus' , nonControlIonEnergy , controlIonEnergy , is_control_qubit);
        
        % 2.3.2 If Acceptable => Increase acceptabilityRate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (depth >= 0) && ((depth - controlRadius) > 0)
            acceptabilityRates(energyLoopID) = acceptabilityRates(energyLoopID) + (1/iterations_per_loop);
        end
    
        % 2.3.3 Progress Probe
        %%%%%%%%%%%%%%%%%%%%%%
        if (mod(currentIteration,50000) == 0)
            test_text_1 = 'In Loop #' + string(energyLoopID) + '  Energy = ' + string(controlIonEnergy) + 'keV  Completed Iterations = ' + string( (currentIteration * 100)/(iterations_per_loop)) + '%';
            disp(test_text_1);
        end
    end

    %{
    % Test
    test_text_2 = 'Loop #' + string(energyLoopIndex) + '  Energy Level = ' + string(controlIonEnergy) + 'keV  =>  Average Depth = ' + string(loop_average) + '  Acceptability Rate = ' + acceptabilityRates(energyLoopIndex)*100 + '%';
    disp(test_text_2);
    %}

end


% 3. ENERGY VS ACCEPTABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on

    % 3.1 INITIALISING FIGURE  
    %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'FontSize',12);

        this_title = "Control Ion Implantation: Ion Energy VS Acceptability Rates";
        title(this_title,'FontWeight','bold','FontSize',18);

        xAxisTitle = "Implantation Energy (keV)";
        yAxisTitle = "Acceptability Rate (%)";

    % 3.2 DISPLAYING FIGURE  
    %%%%%%%%%%%%%%%%%%%%%%%
        xlabel(xAxisTitle, 'color' , 'black' ,'FontSize',14);
        ylabel(yAxisTitle, 'color' , 'black' ,'FontSize',14);
        scatter(energyLevels , acceptabilityRates.*100 , 'red' , 'o' , 'filled' );

hold off

% 4. DEPTH VS ACCEPTABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
    hold on

    % 4.1 INITIALISING FIGURE  
    %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'FontSize',12);

        this_title = "Control Ion Implantation: Average Implantation Depth VS Acceptability Rates";
        title(this_title,'FontWeight','bold','FontSize',18);

        xAxisTitle = "Avarage Implantation Depth (Å)";
        yAxisTitle = "Acceptability Rate (%)";
        
    % 4.2 DISPLAYING FIGURE  
    %%%%%%%%%%%%%%%%%%%%%%%
        xlabel(xAxisTitle, 'color' , 'black' ,'FontSize',14);
        ylabel(yAxisTitle, 'color' , 'black' ,'FontSize',14);
        scatter(averageImplantationDepths , acceptabilityRates.*100 , 'blue' , 'o' , 'filled' );

hold off