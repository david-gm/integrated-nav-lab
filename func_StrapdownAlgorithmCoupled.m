function [ x_e, x_e_spherical,t ] = func_StrapdownAlgorithmCoupled( kinIMUmeassurements,IMUStrapdownConst,constants,...
    kalman,sampleTimeImuForKalman )
%FUNC_STRAPDOWNALGORITHM Summary of this function goes here
%   Detailed explanation goes here

    sampleRate = 50;
    dt = 1/sampleRate;

    t=kinIMUmeassurements(:,1);
    
    
    %rotation (initial-body frame in body frame)
    omega_ib_b=[kinIMUmeassurements(:,3),kinIMUmeassurements(:,2),-kinIMUmeassurements(:,4)];
    
    %specific force(body frame)
    f_b=[kinIMUmeassurements(:,6),kinIMUmeassurements(:,5),-kinIMUmeassurements(:,7)];
    
    [t,omega_ib_b,f_b]=func_Averaging(t,omega_ib_b,f_b,sampleRate);
        
    %correct meassurements with sensor biases
    omega_ib_b = deg2rad(omega_ib_b) - repmat(deg2rad(IMUStrapdownConst.sensorBiases.gyro)/3600,length(t),1);
    f_b = f_b - repmat(IMUStrapdownConst.sensorBiases.acc,length(t),1);

        %preallocation
        Omega_ib_b = zeros(3,3,length(t)+1);
        R_b_L = zeros(3,3,length(t)+1);
        R_b_L_dot = zeros(3,3,length(t)+1);
        Omega_iL_L = zeros(3,3,length(t));
        v_e_L = zeros(3,length(t)+1);
        x_e_spherical = zeros(3,length(t)+1);
        omega_iL_L = zeros(length(t),3);
        f_L = zeros(3,length(t)+1);
        omega_ie_L = zeros(3,length(t));
        g_L = zeros(3,length(t));
        v_dot_e_L = zeros(3,length(t)+1);
        x_dot_e = zeros(3,length(t)+1);
        x_e = zeros(3,length(t)+1);
        
    tk_kalman = 1;    
        
    for tk=1:length(t)
        %% step 1
            % when tk=1, then get R_b_L, x_e_spherical from inital data and set v_e_L=0
            if tk==1
                R_b_L(:,:,1) = getRotationMatrix(deg2rad(IMUStrapdownConst.initialAttitude));
                v_e_L(:,1) = [ 0; 0; 0 ]; 
                
                %initial position
                x_e_spherical(:,1) = IMUStrapdownConst.startingPointRad;
                
                [a_wgs84, b_wgs84] = Projections.getWGS84Axes();
                
                x_e(:,1) = Projections.ellipticToCartesianX(a_wgs84, b_wgs84, x_e_spherical(:,1)')';
                %%%x_e(:,1) = Projections.getCartesianCoords(x_e_spherical(:,1),constants.a,constants.f);
			else 
				% to make sure R_b_L is a rotational matrix: calculate R_b_L from r,p,y;
				R_b_L(:,:,tk) = getRotationMatrix(angles_rpy);
            end
            
            
            % omega_iL_L
            omega_iL_L(tk,:) = get_omega_iL_L(v_e_L(:,tk),x_e_spherical(:,tk),constants);
            
            %==============================================================
            % used values from tk, not from tk+1 cause meassurements start
            % at tk --> correct?
            %==============================================================
            Omega_ib_b(:,:,tk+1) = getSkewSymmetricMatrix(omega_ib_b(tk,:));
            
            Omega_iL_L(:,:,tk) = getSkewSymmetricMatrix(omega_iL_L(tk,:));
                        
            R_b_L_dot(:,:,tk+1) = R_b_L(:,:,tk) * Omega_ib_b(:,:,tk+1) - Omega_iL_L(:,:,tk) * R_b_L(:,:,tk);
                        
            R_b_L_temp = R_b_L(:,:,tk) + dt * R_b_L_dot(:,:,tk+1);
            R_L_b_temp = R_b_L_temp^-1;
			
			% get roll pitch and yaw
			angles_rpy = getRollPitchYaw(R_L_b_temp);
            
        %% step 2
        
            %==============================================================
            % used values from tk, not from tk+1 cause meassurements start
            % at tk --> correct?
            %==============================================================
            f_L(:,tk+1) = R_b_L_temp * f_b(tk,:)';
            
        %% step 3 & 4
        
            [omega_ie_L(:,tk),g_L(:,tk)]=func_getModel_omega_g(x_e_spherical(:,tk),constants);
            Omega_ie_L = getSkewSymmetricMatrix(omega_ie_L(:,tk));
            
            f_coriolis_L = -(Omega_iL_L(:,:,tk) + Omega_ie_L) * v_e_L(:,tk);
            
        %% step 5
        
            %==============================================================
            % should there be g_L(:,tk) OR (like in the paper) g_L(:,tk-1)
            % --> which wouldn't be existent in the first epoch
            %==============================================================
            v_dot_e_L(:,tk+1) = f_L(:,tk+1) + g_L(:,tk) + f_coriolis_L;
            v_e_L(:,tk+1) = v_e_L(:,tk) + dt * v_dot_e_L(:,tk+1);
            
        %% step 6
        
            R_L_e = getR_L_e(x_e_spherical(:,tk));
            
            x_dot_e(:,tk+1) = R_L_e * v_e_L(:,tk+1);
            
            x_e(:,tk+1) = x_e(:,tk) + dt * x_dot_e(:,tk+1);
            
            [a_wgs84, b_wgs84] = Projections.getWGS84Axes();
            
            x_e_spherical(:,tk+1) = Projections.cartesianToEllipticX(a_wgs84,b_wgs84,x_e(:,tk+1)',false);
            %%%x_e_spherical(:,tk+1) = Projections.getSphericalCoords(x_e(:,tk+1),constants.a,constants.f);
            
         %% step 7: kalman filter
         
         % FIXME: kalman jedes mal nicht nur jedes 50ste mal fuer imu,
         % abfrage auch fuer datenluecke dann nur kalman, zaehler tk durch
         % 50 dividieren --> dann erst richtiges t --> imu systemrauschen
         % anpassen --> bei gps wieder absolut also t auf 0 setzen...
         
            % call kalman filter every 50th time --> output x_e and v as
            % input for next epoch of strapdown algorithm
            
            % inputs needed:
            %   1. sampleTimeImuForKalman = 50
            %   2. kalman.in.x_gps       --> should include NaN @ data gaps
            %   3. kalman.init.P_initial
            %   4. kalman.in.precisions
            %   5. kalman.in.dt
            
                % kalman filter
                
                    % define inputs:
                    
                    % for every 50th imu observation, add gps position to
                    % z_k
%                     if mod(tk-sampleTimeImuForKalman/2,sampleTimeImuForKalman)==0 && tk>=sampleTimeImuForKalman/2 
                    if mod(tk,sampleTimeImuForKalman)==0
                        inXGPS = kalman.in.x_gps(:,tk_kalman)'; 

                        tk_kalman = tk_kalman + 1;

                    % for every other epoch put gps observations to NaN
                    else
                        inXGPS = NaN(1,3);
                    end
                    
                    % x_imu: 
                        % is [0 0 0] for 1. epoch and difference for
                        % all other epochs

                    if tk ~= 1
                        kalman.in.x_imu = (x_e(:,tk+1))';% - x_e(:,tk))';                            
                    else
                        kalman.in.x_imu = (x_e(:,tk+1))';%[0 0 0];

                    % x_initial: 
                        % first gps and imu position

                        kalman.init.x_initial = [   kalman.in.x_gps(:,1)
                                                    0;0;0];

                        kalman.out.x_k_hat{1}= kalman.init.x_initial;
                        kalman.out.P_k{1}=kalman.init.P_initial;
                    end

                    % kalman filter
%                     [ kalman.out.x_k_hat{tk+1}, kalman.out.P_k{tk+1} ] = func_kalmanFilteringCoupled( ...
%                         kalman.out.x_k_hat{tk}, ...
%                         kalman.out.P_k{tk}, ...
%                         inXGPS, ...
%                         kalman.in.x_imu, ...
%                         kalman.in.precisions, ...
%                         kalman.in.dt, ...
%                         tk_kalman);
                    [ kalman.out.x_k_hat{tk+1}, kalman.out.P_k{tk+1} ] = func_kalmanFilteringCoupled( ...
                        kalman.out.x_k_hat{tk}, ...
                        kalman.out.P_k{tk}, ...
                        inXGPS, ...
                        kalman.in.x_imu, ...
                        kalman.in.precisions, ...
                        1/50, ...
                        1);
               
                % set output of kalman filter as input of next epoch 
                % for the Strapdown Algorithm    
%                 x_dot_e(:,tk+1) = kalman.out.x_k_hat{tk_kalman+1}(4:6);
%                 x_e(:,tk+1) = kalman.out.x_k_hat{tk_kalman+1}(1:3);
%                 x_e_spherical(:,tk+1) = Projections.getSphericalCoords(x_e(:,tk+1),constants.a,constants.f);    
                        
                x_dot_e(:,tk+1) = kalman.out.x_k_hat{tk+1}(4:6);
                x_e(:,tk+1) = kalman.out.x_k_hat{tk+1}(1:3);
                
                [a_wgs84, b_wgs84] = Projections.getWGS84Axes();
            
                x_e_spherical(:,tk+1) = Projections.cartesianToEllipticX(a_wgs84,b_wgs84,x_e(:,tk+1)',false);
                %%%x_e_spherical(:,tk+1) = Projections.getSphericalCoords(x_e(:,tk+1),constants.a,constants.f);    

    end
        
end

function [t,omega,f] = func_Averaging(t_,omega_,f_,hz_)

    dt=diff(t_);
    hz=1/mean(dt);
    rate=round(hz/hz_);
    
    k=1;
    beginningIndex=1;
    
    for i=1:length(dt)
        if mod(i,rate)==0 
            % change index to i insted of i-1
            % and beginning index to i+1!
            t(k)=(t_(i)+t_(beginningIndex))/2;
            omega(k,:)=mean(omega_(beginningIndex:i,:));
            f(k,:)=mean(f_(beginningIndex:i,:));
            k=k+1;
            beginningIndex=i+1;
        end
    end
        

end

function [R] = getSkewSymmetricMatrix(w)

    R = [ 0    -w(3)  w(2)
          w(3)  0    -w(1)
         -w(2)  w(1)  0   ];
    
end

function [R] = getRotationMatrix(angles)

    % book p.28    
    r=angles(1);
    p=angles(2);
    y=angles(3);

    R_lb = [    cos(p)*cos(y)                         cos(p)*sin(y)                       -sin(p)
                sin(r)*sin(p)*cos(y)-cos(r)*sin(y)    sin(r)*sin(p)*sin(y)+cos(r)*cos(y)  sin(r)*cos(p)
                cos(r)*sin(p)*cos(y)+sin(r)*sin(y)    cos(r)*sin(p)*sin(y)-sin(r)*cos(y)  cos(r)*cos(p)];
    %R_bl       
    R = R_lb^-1;
    
end

function [angles] = getRollPitchYaw(R)

    r=atan2(R(2,3,:),R(3,3,:));
    p=asin(-R(1,3,:));
    y=atan2(R(1,2,:),R(1,1,:));
	
	angles = [r,p,y];

end

function [omega_iL_L] = get_omega_iL_L(v,position,constants)

    phi = position(1);
    h = position(3);
    
    % M and N: curvature radaii
    
        b = constants.a - constants.f * constants.a;
        c = constants.a^2 / b;
        e_line_2 = (constants.a^2 - b^2) / b^2;
                
        V = sqrt(1 + e_line_2 * cos(phi)^2);
        
        M = c/V^3;
        N = c/V;
    
    phi_dot = v(1)/(M+h);
    lambda_dot = v(2)/((N+h) * cos(phi));

    omega_iL_L = [ (lambda_dot + constants.omega_E) * cos(phi)
                   -phi_dot
                   -(lambda_dot + constants.omega_E) * sin(phi)]';
               
end

function [R_L_e] = getR_L_e(position)

   phi = position(1);
   lam = position(2);
    
   n = [-sin(phi)*cos(lam)
        -sin(phi)*sin(lam)
        cos(phi)];
   
   e = [-sin(lam)
         cos(lam)
         0];
    
   d = [-cos(phi)*cos(lam)
        -cos(phi)*sin(lam)
        -sin(phi)];
    
   R_L_e = [n,e,d];

end