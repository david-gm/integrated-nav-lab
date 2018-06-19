function [ omega_ie_l,g_l ] = func_getModel_omega_g( coords,c )
%FUNC_GETMODEL_OMEGA_G Summary of this function goes here
%   Detailed explanation goes here

    phi = coords(1);
    h = coords(3);
    
    %omega_ie_l
    
    omega_ie_l=[c.omega_E*cos(phi)
                0
                -c.omega_E*sin(phi)];
        
    %g_l

        b = c.a-c.f*c.a;
        m = c.omega_E^2 * c.a^2 * b / c.GM_E;
        %formular of somigliana
        sin2phi = sin(phi)^2;
        cos2phi = cos(phi)^2;
        gamma_phi = (c.a * c.gamma_a * cos2phi + b * c.gamma_b * sin2phi) ...
            /(sqrt( c.a^2 * cos2phi + b^2 * sin2phi));
        
        gamma_phi_h = gamma_phi * ( 1 - 2/c.a*(1+c.f+m-2*c.f*sin2phi)*h + 3/c.a^2*h^2);
    
    g_l=[0
        0
        gamma_phi_h];

end

