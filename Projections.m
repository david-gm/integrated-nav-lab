classdef Projections
    %PROJECTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function [A_WGS84, B_WGS84] = getWGS84Axes()
        % getWGS84Axes returns semi-major axis A and semi-minor axis B of
        % WGS84 ellipsoid
        
            f = 1.0/298.257223563;
        
            A_WGS84 = 6378137.0;        % earth semi-major axis (WGS84)
                                        % [m]
                    
            B_WGS84 = (1-f)*A_WGS84;    %earth semi-minor axis (WGS84)
                                        % [m]
            %B_WGS84 = 6356752.3142;     
                                        
        end
                
        function [phi_g, lambda_g, h_g]=cartesianToElliptic(a,b,x,y,z,output_in_deg)
        % cartesianToElliptic
        % INPUTS:	a, b, x, y, z
        %%%%%%%%%%%%%%%%%%%%		description                 datatype
        % a                         major axis                  double
        % b                         semi-major axis             double
        % x, y and z                cartesian coordinates       double
        % output_in_deg             returns degrees if true     bool
        % OUTPUTS:  phi_g, lambda_h, h_g in radians/degrees
        %%%%%%%%%%%%%%%%%%%%		description                 datatype
        % phi_g, lambda_h, h_g		ellipsoidical coordinates	double

            p=sqrt(x.^2+y.^2);
            e_2=(a^2-b^2)/a^2;
      
            lambda_g = atan2(y,x);
            phi_g1 = atan2(z,p*(1-e_2));
            
            % 0.001 sec
            precision = deg2rad(10^-6/3600);
            
            for i=1:length(x)
                diffPhi_g = 1;
                while (diffPhi_g > precision)
                    N(i)=a^2/sqrt(a^2*cos(phi_g1(i))^2 + ...
                        b^2*sin(phi_g1(i))^2);
                    h_g(i) = p(i)/cos(phi_g1(i))-N(i);
                    phi_g2(i) = atan2(z(i), p(i) * ...
                        (1 - e_2 * N(i)/(N(i)+h_g(i))));

                    diffPhi_g = abs(phi_g1(i)-phi_g2(i));
                    phi_g1(i) = phi_g2(i);
                end
                phi_g(i) = phi_g1(i);
            end
            
            [n,m] = size(x);
            if n>m
                phi_g = phi_g';
                h_g = h_g';
            end
            
            if output_in_deg
                phi_g = rad2deg(phi_g);
                lambda_g = rad2deg(lambda_g);
            end
        end

        function [ell_pos]=cartesianToEllipticX(a,b,X,output_in_deg)
        % cartesianToEllipticX
        % INPUTS:   a, b, X
        %%%%%%%%%%%%%%%%%%%%        description                 datatype
        % a                         major axis                  double
        % b                         semi-major axis             double
        % X                         cartesian coordinates       matrix[double]
        % output_in_deg             returns degrees if true     bool
        % OUTPUTS:  ell_pos in [radians/degrees, m]
        %%%%%%%%%%%%%%%%%%%%        description                 datatype
        % ell_pos                   ellipsoidical coordinates   matix[double]

            size_x = size(X);
            if size_x(1) == 3
                x = X(1,:);
                y = X(2,:);
                z = X(3,:);
            elseif size_x(2) == 3
                x = X(:,1)';
                y = X(:,2)';
                z = X(:,3)';
            else
                error('Provided matrix X is not of size(x,3) or size(3,x)')
            end
            p=sqrt(x.^2+y.^2);
            e_2=(a^2-b^2)/a^2;
      
            lambda_g = atan2(y,x);
            phi_g1 = atan2(z,p*(1-e_2));
            
            % 0.001 sec
            precision = deg2rad(10^-6/3600);
            
            for i=1:length(x)
                diffPhi_g = 1;
                while (diffPhi_g > precision)
                    N(i)=a^2/sqrt(a^2*cos(phi_g1(i))^2 + ...
                        b^2*sin(phi_g1(i))^2);
                    h_g(i) = p(i)/cos(phi_g1(i))-N(i);
                    phi_g2(i) = atan2(z(i), p(i) * ...
                        (1 - e_2 * N(i)/(N(i)+h_g(i))));

                    diffPhi_g = abs(phi_g1(i)-phi_g2(i));
                    phi_g1(i) = phi_g2(i);
                end
                phi_g(i) = phi_g1(i);
            end
            
            
            if output_in_deg
                phi_g = rad2deg(phi_g);
                lambda_g = rad2deg(lambda_g);
            end
                        
            if size_x(1) == 3
                ell_pos = [phi_g; lambda_g; h_g];
            else
                ell_pos = [phi_g', lambda_g', h_g'];
            end  

        end
        
        function [x,y,z] = ellipticToCartesian(a,b,phi,lam,h)
        % ellipticToCartesian
        % INPUTS:	a, b, phi, lam, h
        %%%%%%%%%%%%%%%%%%%%		description                   datatype
        % a                         major axis                    double
        % b                         semi-major axis               double
        % phi, lam and h            elliptic coordinates(radians) double
        % OUTPUTS:  x, y, z
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % x, y, z                   cartesian coordinates   double
            
            N=a^2./sqrt(a^2*cos(phi).^2 + b^2*sin(phi).^2);
            
            x = (N+h) .* cos(phi) .* cos(lam);
            y = (N+h) .* cos(phi) .* sin(lam);
            z = (b^2/a^2 * N + h) .* sin(phi);
        end
        
        function [X] = ellipticToCartesianX(a,b,ellip_coords)
        % ellipticToCartesianX
        % INPUTS:	a, b, ellip_coords
        %%%%%%%%%%%%%%%%%%%%		description                         datatype
        % a                         major axis                          double
        % b                         semi-major axis                     double
        % ellip_coords              elliptic coordinates(radians, m)    matrix: [n x [phi, lam, h]]     double
        % OUTPUTS:  X
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % X: [x, y, z]              cartesian coordinates   double
            phi = ellip_coords(:,1);
            lam = ellip_coords(:,2);
            h = ellip_coords(:,3);
            
            N=a^2./sqrt(a^2*cos(phi).^2 + b^2*sin(phi).^2);
            
            x = (N+h) .* cos(phi) .* cos(lam);
            y = (N+h) .* cos(phi) .* sin(lam);
            z = (b^2/a^2 * N + h) .* sin(phi);
            
            X = [x,y,z];
        end
        
        function [x,y,z] = ellipticToCartesianAsInICD(a,e2,phi,lam,h)
        % ellipticToCartesian
        % INPUTS:	a, b, phi, lam, h
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % a                         major axis              double
        % e2                        first numerical eccentricity  double
        % phi, lam and h            elliptic coordinates    double
        % OUTPUTS:  x, y, z
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % x, y, z                   cartesian coordinates   double
            
            N = a ./ sqrt(1 - e2 * sin(phi).^2);
            
            x = (N+h) .* cos(phi) .* cos(lam);
            y = (N+h) .* cos(phi) .* sin(lam);
            z = ((1 - e2) * N + h) .* sin(phi);
        end
        
        function [a,z] = cartesianToLLF_AZ(dX,dY,dZ,phi,lam)
        % cartesianToLLF_AZ
        % INPUTS:	dX,dY,dZ,phi,lam
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % dX                        difference vector [m]   double
        % dY
        % dZ
        % phi                       location of user  [rad] double
        % lam                       
        % OUTPUTS:  a, z
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % a                         azimuth [rad]           double
        % z                         zenith distance [rad]   double
        
            % LL north
            n = -dX .* sin(phi) .* cos(lam) - ...
                 dY .* sin(phi) .* sin(lam) + ...
                 dZ .* cos(phi);
            % LL east
            e = -dX .*sin(lam) + dY .*cos(lam);
            % LL up
            u =  dX .* cos(phi) .* cos(lam) + ...
                 dY .* cos(phi) .* sin(lam) + ...
                 dZ .* sin(phi);
             
            % azimuth
            a = atan2(e,n);
            % zenith distance
            z = acos(u./sqrt(n.^2 + e.^2 + u.^2));          
            
        end
        
        function pos = cartesianToLLF(dX,dY,dZ,phi,lam)
        % cartesianToLLF
        % INPUTS:	dX,dY,dZ,phi,lam
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % dX                        difference vector [m]   double
        % dY
        % dZ
        % phi                       location of user  [rad] double
        % lam                       
        % OUTPUTS:  pos
        %%%%%%%%%%%%%%%%%%%%		description             datatype
        % pos                       [x[m],y[m],z[m]]        vector(x,3) or vector(3,x)
        %                           in LLF cartesian coordinates
            % LL north
            n = -dX .* sin(phi) .* cos(lam) - ...
                 dY .* sin(phi) .* sin(lam) + ...
                 dZ .* cos(phi);
            % LL east
            e = -dX .*sin(lam) + dY .*cos(lam);
            % LL up
            u =  dX .* cos(phi) .* cos(lam) + ...
                 dY .* cos(phi) .* sin(lam) + ...
                 dZ .* sin(phi);
            
            size_x = size(dX);
                         
            if size_x(1) > size_x(2)
                pos = [n,e,u];
            else
                pos = [n',e',u'];
            end

        end
        
        function pos = LLToCartesian(dx,dy,dz,phi,lam)
            % LLToCartesian
            % INPUTS:	dy,dy,dz,phi,lam
            %%%%%%%%%%%%%%%%%%%%		description             datatype
            % dx                        difference vector [m]   double
            % dy                        in LLF
            % dz
            % phi                       location of user  [rad] double
            % lam                       
            % OUTPUTS:  pos
            %%%%%%%%%%%%%%%%%%%%		description             datatype
            % pos                       [X[m],Y[m],Z[m]]        vector(i,3)
            %                           in global cartesian coordinates
            
            R = [ -sin(phi)*cos(lam)    -sin(lam)   cos(phi)*cos(lam)
                  -sin(phi)*sin(lam)     cos(lam)   cos(phi)*sin(lam)
                   cos(phi)              0          sin(phi)
                ];
            
            pos = zeros(length(dx), 3);
            
            for i=1:length(dx)
                x = [dx(i), dy(i), dz(i)]';
                pos(i,:) = R * x;
            end
        end
        
        function pos = transverseMercatorToEllipsoidalCoords(a,b,x,y,lam0)
            % TRANSVERSEMERCATORTOELLIPSOIDALCOORDS
                        
            % define auxiliary variables
            
            n = (a-b)/(a+b);
            a_ = (a+b)/2 * (1 + 1/4*n^2 + 1/64*n^4);
            b_ = 3/2*n - 27/32*n^3 + 269/512*n^5;
            c_ = 21/16*n^2 - 55/32*n^4;
            d_ = 151/96*n^3 - 417/128*n^5;
            e_ = 1097/512*n^4;
            
            %%%y_ = y/a_;
            x_ = x/a_;
            
            % get fingerprint latitude
            %%%phi_f = y_ + b_*sin(2*y_) + c_*sin(4*y_) + d_*sin(6*y_) + ...
            %%%   e_*sin(8*y_);
            phi_f = x_ + b_*sin(2*x_) + c_*sin(4*x_) + d_*sin(6*x_) + ...
                e_*sin(8*x_);
                
            
%             % arc length of meridian
%             B = a_ * (phi_f + b_*sin(2*phi_f) + c_*sin(4*phi_f) + ...
%                 d_*sin(6*phi_f) + e_*sin(8*phi_f));

            % radius curvature in prime vertical
            N = a^2 / sqrt(a^2*cos(phi_f)^2 + b^2*sin(phi_f)^2);
            
            % auxilary quantities
            eta = cos(phi_f)/b * sqrt(a^2 - b^2);
            t = tan(phi_f);
            
            % ellipsoidal coords
            phi = phi_f + t/(2*N^2)*(-1-eta^2)*y^2 + ...
                t/(24*N^4) * (5 + 3*t^2 + 6*eta^2 - 6*t^2*eta^2 - ...
                    3*eta^4 - 9*t^2*eta^4) * y^4 + ...
                t/(720*N^6) * (-61 - 90*t^2 - 45*t^4 - 107*eta^2 + ...
                    162*t^2*eta^2 + 45*t^4+eta^2) * y^6 + ...
                t/(40320*N^8) * (1385 + 3633*t^2 + 4095*t^4 + 1575*t^6)...
                * y^8;
                                 
            lambda = lam0 + 1/(N*cos(phi_f))*y + 1/(6*N^3*cos(phi_f)) * ...
                        (-1 - 2*t^2 - eta^2) * y^3 + ...
                     1/(120*N^5*cos(phi_f)) * (5 + 28*t^2 + 24*t^4 + ...
                        6*eta^2 + 8*t^2*eta^2) * y^5 + ...
                     1/(5040*N^7*cos(phi_f)) * (-61 -662*t^2 - ...
                     1320*t^4 - 720*t^6) * y^7;

             pos = [phi, lambda];
            
        end
        
        function pos = ellipsoidalToTransverseMercatorCoords(a,b,phi,lam,lam0)
            % ELLIPSOIDALTOTRANSVERSEMERCATORCOORDS
            
            % define auxiliary variables
            
            n = (a-b)/(a+b);
            a_ = (a+b)/2 * (1 + 1/4*n^2 + 1/64*n^4);
            b_ = - 3/2*n - 9/16*n^3 - 3/32*n^5;
            c_ = 15/16*n^2 - 15/32*n^4;
            d_ = -35/48*n^3 + 105/256*n^5;
            e_ = 315/512*n^4;
            
            % arc length of meridian
            B = a_*(phi + b_*sin(2*phi) + c_*sin(4*phi) + d_*sin(6*phi) ...
                + e_*sin(8*phi));
            
            % radius curvature in prime vertical
            N = a^2 / sqrt(a^2*cos(phi)^2 + b^2*sin(phi)^2);
            
            % auxilary quantities
            eta = cos(phi)/b * sqrt(a^2 - b^2);
            t = tan(phi);
            
            % longitude difference
            l = lam - lam0;
            
            % plane coords
            y = B + t/2*N*cos(phi)^2*l^2 + ...
                t/24*N*cos(phi)^4*(5-t^2+9*eta^2+4*eta^4)*l^4 + ...
                t/720*N*cos(phi)^6*...
                    (61-58*t^2+t^4+270*eta^2-330*t^2*eta^2)*l^6 + ...
                t/40320*N*cos(phi)^8*...
                    (1385-3111*t^2+543*t^4-t^6)*l^8;
                
            x = N*cos(phi)*l + 1/6*N*cos(phi)^3*(1-t^2+eta^2)*l^3 + ...
                1/120*N*cos(phi)^5*...
                    (5-18*t^2+t^4+14*eta^2-58*t^2*eta^2)*l^5 + ...
                1/5040*N*cos(phi)^7*(61-479*t^2+179*t^4-t^6)*l^7;
            
            pos = [x,y];
        end
        
        function pos = gaussKruegerEPSG31256ToWGS84Cart(lam0, x, y, H, N, offset)
            % GAUSSKRUEGEREPSG31256TOWGS84
            
            % Bessel Parameter a,b:
            A_BESSEL = 6377397.15508;   % earth semi-major axis (Bessel)
                                        % [m]
            
            B_BESSEL = 6356078.9629;    %earth semi-minor axis (Bessel)
                                        % [m]
            
            % apply offset to x coordinate
            x = x+offset;
                                        
            % planar coords to phi, lam, h on Bessel Ellipsoid
            phiLam = Projections.transverseMercatorToEllipsoidalCoords(...
                A_BESSEL,B_BESSEL,x,y,lam0);
            
            % apply geoid undulation N
            h = N+H;
            
            % phi, lam, h to X,Y,Z
            [xBessel,yBessel,zBessel] = Projections.ellipticToCartesian(...
                A_BESSEL, B_BESSEL, phiLam(1), phiLam(2), h);

            % transformation Bessel -> WGS84 Ellipsoid
            
                % definition of parameters (WGS->Ellipsoid)
                % cx = -577.33; cy = -90.13; cz = -463.92; [m] -> t
                % d_1 = 15.864; d_2 = 4.537; d_3 = 16.358; [cc]-> alpha
                % d = -2.42*10^-6; [ppm] -> mu

				% script of reference systems
				transParWGS84_Bessel.t = [-577.33; -90.13; -463.92];
                
                % al from script of reference systems
                
                transParWGS84_Bessel.al = [15.8537; 4.5500; 16.3489] ...
                    / 10000 / 400 * 2*pi;
                
%                 transParWGS84_Bessel.al = [15.864; 4.537; 16.358] ...
%                     / 10000 / 400 * 2*pi;

				% script of reference systems
                transParWGS84_Bessel.mu = 1-2.42*10^-6;

                pos = Transformations.revHelmert7(...
                    [xBessel,yBessel,zBessel]', ...
                    transParWGS84_Bessel.al, ...
                    transParWGS84_Bessel.t, ...
                    transParWGS84_Bessel.mu);
        end
        
        function pos = wgs84CartToGaussKruegerEPSG31256(lam0, x, y, z, N, offset)
            % WGS84CARTTOGAUSSKRUEGEREPSG31256
            
            % Bessel Parameter a,b:
            A_BESSEL = 6377397.15508;   % earth semi-major axis (Bessel)
                                        % [m]
            
            B_BESSEL = 6356078.9629;    %earth semi-minor axis (Bessel)
                                        % [m]
            
  
            wgs84Cart_ = [x,y,z]';

            % Transformation WGS84 -> Bessel
                transParWGS84_Bessel.t = [-577.33; -90.13; -463.92];              
                
                transParWGS84_Bessel.al = [15.864; 4.537; 16.358] ...
                    / 10000 / 400 * 2*pi;

                transParWGS84_Bessel.mu = 1-2.42*10^-6;

                besselCart = Transformations.helmert7(...
                    wgs84Cart_, ...
                    transParWGS84_Bessel.al, ...
                    transParWGS84_Bessel.t, ...
                    transParWGS84_Bessel.mu);

            % Cart Bessel -> phi, lam, h
            [phi,lam,h] = Projections.cartesianToElliptic(...
                A_BESSEL, B_BESSEL, ...
                besselCart(1), besselCart(2), besselCart(3));
            besselElliptic = [phi,lam,h];

            % Bessel phi, lam, h -> Projection to GK M34 East

            pos = Projections.ellipsoidalToTransverseMercatorCoords(...
                A_BESSEL,B_BESSEL,...
                besselElliptic(1),besselElliptic(2),lam0);

            % apply negative offset to x
            pos(2) = pos(2)-offset;

            H = h-N;

            pos = [pos, H];
            
        end
    end
    
end

