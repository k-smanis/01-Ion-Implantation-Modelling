% This function returns the trajectorial change that the Phosphorus Ion undergoes during its implantation in Silicon
% Depending on the input, the output is the change in the x or the y or the z axis

function y = path_along_axis(axis , ionType , nonControlIonEnergy, controlIonEnergy , is_control_qubit)

% 1. TESTING FOR CONTROL QUBIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (is_control_qubit)
    ionType = 'P';
    ionEnergy = controlIonEnergy;
else
    ionEnergy = nonControlIonEnergy;
end
    
% 2. SETTING POLYNOMIAL REGRESSION COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (string(ionType) == 'P')
        if (axis == 'x')
            s1 = (-1.8956805899)*(10^(-13)); s2 = (3.3923618047)*(10^(-10)); s3 = (-2.4645391336)*(10^(-7)); s4 = (9.5394043242)*(10^(-5)); s5 = (-2.3643081247)*(10^(-2)); s6 = (6.0991747094)*(10^(0)); straggling_const = (2.2054301983)*(10^(1));
            r1 = (4.0745059722)*(10^(-14)); r2 = (-3.3601039663)*(10^(-11)); r3 = (-7.3049769110)*(10^(-9)); r4 = (1.1186615874)*(10^(-5)); r5 = (-4.2105628811)*(10^(-3)); r6 = (1.3398417218)*(10^(1)); range_const = (3.2060314440)*(10^(1));
        elseif (axis == 'y') || (axis == 'z')
            s1 = (-2.3868298905)*(10^(-13)); s2 = (3.7852777809)*(10^(-10)); s3 = (-2.3439710828)*(10^(-7)); s4 = (7.3120751592)*(10^(-5)); s5 = (-1.3926123432)*(10^(-2)); s6 = (4.1931621904)*(10^(0)); straggling_const = (1.7594037870)*(10^(1));
            r1 = 0; r2 = 0; r3 = 0; r4 = 0; r5 = 0; r6 = 0; range_const = 0;
        end
    elseif (string(ionType) == 'Sb')
        if (axis == 'x')
            s1 = (-3.2690443147)*(10^(-13)); s2 = (5.2448917764)*(10^(-10)); s3 = (-3.2425963991)*(10^(-7)); s4 = (9.7453022292)*(10^(-5)); s5 = (-1.5078499426)*(10^(-2)); s6 = (2.0053312158)*(10^(0)); straggling_const = (1.8591837191)*(10^(1));
            r1 = (-8.8565105611)*(10^(-13)); r2 = (1.4252712448)*(10^(-9)); r3 = (-8.8606740568)*(10^(-7)); r4 = (2.6883211211)*(10^(-4)); r5 = (-4.1682800162)*(10^(-2)); r6 = (7.1589482290)*(10^(0)); range_const = (4.7729568585)*(10^(1));
        elseif (axis == 'y') || (axis == 'z')
            s1 = (-2.2159096507)*(10^(-13)); s2 = (3.5662233870)*(10^(-10)); s3 = (-2.2285015088)*(10^(-7)); s4 = (6.8543650969)*(10^(-5)); s5 = (-1.1061033432)*(10^(-2)); s6 = (1.6015728512)*(10^(0)); straggling_const = (1.3172298862)*(10^(1));
            r1 = 0; r2 = 0; r3 = 0; r4 = 0; r5 = 0; r6 = 0; range_const = 0;
        end
    elseif (string(ionType) == 'Bi')
        if (axis == 'x')
            s1 = (-2.7982178729)*(10^(-13)); s2 = (4.4777755764)*(10^(-10)); s3 = (-2.7779161108)*(10^(-7)); s4 = (8.4628926911)*(10^(-5)); s5 = (-1.3465236936)*(10^(-2)); s6 = (1.4917409794)*(10^(0)); straggling_const = (1.6479714742)*(10^(1));
            r1 = (-9.6413455226)*(10^(-13)); r2 = (1.5479457319)*(10^(-9)); r3 = (-9.6481992870)*(10^(-7)); r4 = (2.9614565783)*(10^(-4)); r5 = (-4.7774657090)*(10^(-2)); r6 = (6.4514781640)*(10^(0)); range_const = (5.6604602995)*(10^(1));
        elseif (axis == 'y') || (axis == 'z')
            s1 = (-1.8422206816)*(10^(-13)); s2 = (2.9529194422)*(10^(-10)); s3 = (-1.8366687079)*(10^(-7)); s4 = (5.6311745156)*(10^(-5)); s5 = (-9.1486015524)*(10^(-3)); s6 = (1.1404020764)*(10^(0)); straggling_const = (1.1444802024)*(10^(1));
            r1 = 0; r2 = 0; r3 = 0; r4 = 0; r5 = 0; r6 = 0; range_const = 0;
        end
    end


% 3. RANGE CALCULATION
%%%%%%%%%%%%%%%%%%%%%%
    range_function = @(x) r1 * x.^6 + r2 * x.^5 + r3 * x.^4 + r4 * x.^3 + r5 * x.^2 + r6 * x + range_const;
    range = range_function(ionEnergy);

% 4. STRAGGLING CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    straggling_function = @(x) s1 * x.^6 + s2 * x.^5 + s3 * x.^4 + s4 * x.^3 + s5 * x.^2 + s6 * x + straggling_const;
    straggling = straggling_function(ionEnergy);


% 5. Monte Carlo Validity Test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    if (axis == 'x')
        distributionTitle = "Longitudinal Distribution";
        xAxisTitle = "Longitudinal Depth (Å)";
    elseif (axis == 'y') || (axis == 'z')
        distributionTitle = "Lateral Distribution";
        xAxisTitle = "Lateral Spread (Å)";
    end
    DATA = normrnd(range,straggling , 5000 , 1);
    figure
        hold on;
            title(distributionTitle);
            histogram(DATA,200);
            histfit(DATA);
            ImplantationDataDistribution = fitdist(DATA,'Normal');
            xlabel(xAxisTitle);

            annotationText1 = "From Plotted Data:" + newline + "Mean = " + string(mean(ImplantationDataDistribution)) + 'Å' + newline + "Straggling = " + string(std(ImplantationDataDistribution)) + 'Å';
            annotationText2 = "From Polynomial Regression:" + newline + "Mean = " + string(range) + 'Å' + newline + "Straggling = " + string(straggling) + 'Å';
            annotation('textbox',[0.166032210834556 0.806451612903226 0.179355783308932 0.0937019969278033],'Color','red','VerticalAlignment','middle','String',annotationText1,'HorizontalAlignment','center','FontWeight','bold','FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8]);
            annotation('textbox',[0.703513909224016 0.806451612903226 0.179355783308932 0.0937019969278033],'Color','blue','VerticalAlignment','middle','String',annotationText2,'HorizontalAlignment','center','FontWeight','bold','FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8]);

        hold off;
%}


% 6. Polynomial Regression Validity Test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
text = string(ionType) + ', ' + string(ionEnergy)+ 'keV, ' + string(axis) + '-axis   =>   Range=' + string(range) + ' Å    Straggling=' + string(straggling) +' Å';
disp(text);
%}

% 7. RETURNING PATH ALONG SPECIFIED AXIS
% normrnd(μ,σ) returns a random number under the normal probabily 
% distribution with mean = μ   standard deviation = σ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = normrnd(range,straggling);

end