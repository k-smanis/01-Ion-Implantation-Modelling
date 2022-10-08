function print_tragectories(ionType , nonControlIonEnergy , controlIonEnergy , currentSpacing , controlRadius , nonControlRadius , ImpactLocations , DestinationLocations)
figure;
hold on;

    % 1. Initialise Figure Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axis equal;
        implantation_title = 'Implantation Parameters:  Control Ion(Phosphorus,' + string(controlIonEnergy) + 'keV)  Non-Control Ion(' + string(ionType) + ','+ string(nonControlIonEnergy) + 'keV)' + '  Spacing(' + string(currentSpacing) + 'Å)';
        title(implantation_title);
        xlabel('X - Longitudinal Range (Å)');
        ylabel('Y - Lateral Range (Å)');
        zlabel('Z - Lateral Range (Å)');
        view(3);
        grid on;

    % 2. Silicon Surface Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        patch([0,0,0,0],currentSpacing*[-2,6,6,-2],currentSpacing*[-2,-2,6,6],'black', 'FaceAlpha' , 0.5);
   
    % 3. Array Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1 : 9
            % 3.1 Ion Paths
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            scatter3(ImpactLocations(i,1),ImpactLocations(i,2),ImpactLocations(i,3), 50, 'o', 'filled', 'red');
            scatter3(DestinationLocations(i,1),DestinationLocations(i,2),DestinationLocations(i,3), 10,'d', 'filled', 'green');
            line([ImpactLocations(i,1), DestinationLocations(i,1);], [ImpactLocations(i,2), DestinationLocations(i,2);],[ImpactLocations(i,3),DestinationLocations(i,3)] ,'LineWidth', 2);
    
            % 3.2 Wave Functions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (i==5)
                [X,Y,Z] = sphere;
                X1 = X * controlRadius;
                Y1 = Y * controlRadius;
                Z1 = Z * controlRadius;
                surf(X1+DestinationLocations(i,1),Y1+DestinationLocations(i,2),Z1+DestinationLocations(i,3),'FaceAlpha', 0.3, 'EdgeColor','none', 'FaceColor' , 'flat');
            else
                [X,Y,Z] = sphere;
                X1 = X * nonControlRadius;
                Y1 = Y * nonControlRadius;
                Z1 = Z * nonControlRadius;
                surf(X1+DestinationLocations(i,1),Y1+DestinationLocations(i,2),Z1+DestinationLocations(i,3),'FaceAlpha', 0.3, 'EdgeColor','none', 'FaceColor' , 'flat');
            end
        end

hold off;
end