function y = spherical_overlap(DestinationLocations , ion1 , ion2 , nonControlRadius , controlRadius , ionPairID)

% 1. INITIALISING SPHERE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%
    s1_x = DestinationLocations(ion1, 1);
    s1_y = DestinationLocations(ion1, 2);
    s1_z = DestinationLocations(ion1, 3);

    if (ion1 == 5) || (ion2 ==5)
        r1 = controlRadius;
        isControlOverlap = true;
    else
        r1 = nonControlRadius;
        isControlOverlap = false;
    end

% 2. INITIALISING SPHERE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%
    s2_x = DestinationLocations(ion2, 1);
    s2_y = DestinationLocations(ion2, 2);
    s2_z = DestinationLocations(ion2, 3);

    r2 = nonControlRadius;

% 3. SETTING FUNCTION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distance = sqrt( ( s1_x - s2_x )^2 + ( s1_y - s2_y )^2 +( s1_z - s2_z )^2 );

    if ( (distance-r2) > r1)
        tooFarToOverlap = true;
    else
        tooFarToOverlap = false;
    
        if (~isControlOverlap) || (distance == 0)
            x_intercept = distance / 2;
        elseif(isControlOverlap)
            x_intercept = (distance^2 + r1^2 - r2^2) / (2*distance);
            k = distance - x_intercept;
        end
    end

% 4. CALCULATING SPHERICAL OVERLAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if(tooFarToOverlap)
        % Case 1: No Overlap (Too far for it)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        overlap = 0;
    elseif (isControlOverlap)
        
        % Case 2: The nonControl is completely contained in the Control
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((distance + r2) < r1)
            overlap = (4*pi/3) * r2^3;
            mostlyOverlapping = false;
            mostlyNotOverlapping = false;
            completelyOverlapping = true;

        % Case 3: The nonControl is "mostly out of" the Control
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif (x_intercept <= distance)
            mostlyOverlapping = false;
            mostlyNotOverlapping = true;
            completelyOverlapping = false;
            % Control Cap
            fnctn1 = @(phi,theta,rho) rho.^2 .* sin(phi);
            phi_min = 0; phi_max = pi; theta_min = -pi/2; theta_max = pi/2; 
            rho_min = 0; rho_max1 = @(phi,theta) (-x_intercept*sin(phi).*cos(theta)) + sqrt((x_intercept)^2 * (((sin(phi)).^2) .* ((cos(theta)).^2) - 1) + r1^2);
            % s1Cap Volume
            s1CapVolume_1 = integral3(fnctn1,phi_min,phi_max,theta_min,theta_max,rho_min,rho_max1,'Method','tiled');

            % nonControl Cap
            fnctn2 = @(phi,theta,rho) rho.^2 .* sin(phi);
            phi_min = 0; phi_max = pi; theta_min = pi/2; theta_max = 3*(pi/2);
            rho_min = 0; rho_max2 = @(phi,theta) (k*sin(phi).*cos(theta)) + sqrt( k^2 * (((sin(phi)).^2) .* ((cos(theta)).^2 ) - 1) + r2^2);
            % s2Cap Volume
            s2CapVolume_1 = integral3(fnctn2,phi_min,phi_max,theta_min,theta_max,rho_min,rho_max2,'Method','tiled');
    
            % Total Overlap
            overlap = s1CapVolume_1 + s2CapVolume_1;

        % Case 4: The nonControl is "mostly in" the Control
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif (x_intercept > distance)
            mostlyOverlapping = true;
            mostlyNotOverlapping = false;
            completelyOverlapping = false;
            % Control Cap
            fnctn3 = @(phi,theta,rho) rho.^2 .* sin(phi);
            theta_min = -pi/2; theta_max = pi/2; phi_min = 0; phi_max = pi; 
            rho_min = 0; rho_max3 = @(phi,theta) (-x_intercept*sin(phi).*cos(theta)) + sqrt((x_intercept)^2 * (((sin(phi)).^2) .* ((cos(theta)).^2) - 1) + r1^2);
            % s1Cap Volume
            s1CapVolume_2 = integral3(fnctn3,phi_min,phi_max,theta_min,theta_max,rho_min,rho_max3,'Method','tiled');

            % nonControl Cap
            fnctn4 = @(phi,theta,rho) rho.^2 .* sin(phi);
            theta_min = -pi/2; theta_max = pi/2; phi_min = 0; phi_max = pi;
            rho_min = 0; rho_max4 = @(phi,theta) (k*sin(phi).*cos(theta)) + sqrt( k^2 * (((sin(phi)).^2) .* ((cos(theta)).^2) - 1) + r2^2);
            % _negative_s2Cap Volume
            negative_s2CapVolume_2 = ((4*pi/3) * r2^3)-(integral3(fnctn4,phi_min,phi_max,theta_min,theta_max,rho_min,rho_max4,'Method','tiled'));
    
            % Total Overlap
            overlap = s1CapVolume_2 + negative_s2CapVolume_2;
        end

    elseif (~isControlOverlap)
        mostlyOverlapping = false;
        mostlyNotOverlapping = false;
        completelyOverlapping = false;
        
        % Case 1: Spheres completely overlap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (distance == 0)
            overlap = (4*pi/3) * r2^3;
            
        % Case 2: Spheres partially overlap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            fnctn5 = @(phi,theta,rho) rho.^2 .* sin(phi);
            theta_min = -pi/2; theta_max = pi/2; phi_min = 0; phi_max = pi;
            rho_min = 0; rho_max5 = @(phi,theta) (-x_intercept*sin(phi).*cos(theta)) + sqrt((x_intercept)^2 * (((sin(phi)).^2) .* ((cos(theta)).^2) - 1) + r2^2);
    
            % Total Overlap (2 identical caps)
            overlap = 2 * integral3(fnctn5,phi_min,phi_max,theta_min,theta_max,rho_min,rho_max5,'Method','tiled');
        end

    end

% 5. VERIFICATION TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if (ionPairID > 0)
    text0 = 'Possible Spherical Overlap #' + string(ionPairID) + ': (Ion#' + string(ion1) + ') with (Ion#' + string(ion2) + ') :';
    text1 = 'Sphere 1('+string(s1_x)+','+string(s1_y)+','+string(s1_z)+').....Sphere 2('+string(s2_x)+','+string(s2_y)+','+string(s2_z);
    text2 = 'Overlap Volume = ' + string(overlap) + ' Å^3' + newline;
    if (tooFarToOverlap)
        text3 = 'Distance='+string(distance) + '......Too far to overlap';
    elseif (mostlyOverlapping)
        text3 = 'Distance='+string(distance) + '......Mostly overlapping';
    elseif (mostlyNotOverlapping)
        text3 = 'Distance='+string(distance) + '......Mostly not overlapping';
    elseif (~isControlOverlap)
        text3 = 'Distance='+string(distance) + '......Overlapping';
    elseif (completelyOverlapping)
        text3 = 'Distance='+string(distance) + '......Completely Overlapping';
    end

    disp( text0 );disp( text1 ); disp( text2 ); disp( text3 );
end
%}

% 6. RETURNING OVERLAP
%%%%%%%%%%%%%%%%%%%%%%
y = overlap;
end