function [chi_d, U_d] = los(WP, p, p0)
    % Parameters
    L_pp = 304.8; % [m]
    U_d = 6.63; % [m/s]
    Delta = 2*L_pp;
    R = L_pp;

    persistent counter;
    persistent i; % Index of waypoint to go to
    persistent p_k;
    persistent p_kp1;
    
    if isempty(counter)
        counter = 1;
        
        i = 1;
        p_k = p0;
        p_kp1 = WP(:,i);
    end 
    
    if counter == 2
        p_k = p;
        p_kp1 = WP(:,i);
    elseif norm(p - p_kp1) < R
        if i == size(WP,2)
            counter = 2;
        else
            i = i + 1;
        end
        p_k = WP(:,i-1);
        p_kp1 = WP(:,i);
    end
    
    alpha_k = atan2(p_kp1(2) - p_k(2), p_kp1(1) - p_k(1));
    R = [cos(alpha_k), -sin(alpha_k);
        sin(alpha_k), cos(alpha_k)];
    epsilon = R'*(p - p_k);
    e = epsilon(2);
    
    chi_p = alpha_k;
    chi_r = atan2(-e , Delta);
    %rad2deg(chi_p)
    %rad2deg(chi_r)
    chi_d = chi_p + chi_r;
    
end

