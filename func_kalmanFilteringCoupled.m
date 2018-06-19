function [ x_k_hat, P_k ] = func_kalmanFilteringCoupled( x_k_hat, P_k, x_gps, x_imu, prec, dt, t )
%FUNC_KALMANFILTERING Summary of this function goes here
%   Detailed explanation goes here  

        Phi_k=[1 0 0 dt 0 0
                  0 1 0 0 dt 0
                  0 0 1 0 0 dt
                  0 0 0 1 0 0
                  0 0 0 0 1 0
                  0 0 0 0 0 1];

        %sx=(0.5*1.4*10^(-3)*t^2+1/6*3.8*10^(-8)*t^3);
        %sx = (0.5*1.4*10^(-3)*t^2+1/6*3.8*10^(-8)*t^3)*10^5;
        
        if ~isnan(x_gps)
%             t = 20;
            t = 2*10^4;
            sx = (0.5*1.4*10^(-3)*t^2+1/6*3.8*10^(-8)*t^3);
            
            H_k=[eye(3),zeros(3);eye(3),zeros(3)];
            R_k=diag([prec.gps^2;prec.gps^2;prec.gps^2;sx^2;sx^2;sx^2]);
            
            % Beobachtungsvektor
            z_k=[ x_gps'
                  x_imu'];
              
              
        else
            sx = (0.5*1.4*10^(-3)*t^2+1/6*3.8*10^(-8)*t^3);
            
            H_k=[eye(3),zeros(3)];
            R_k=diag([sx^2;sx^2;sx^2]);
            
            
            % Beobachtungsvektor
            z_k=x_imu';
        end

        
        N_k=[1/2*dt^2 0 0
               0 1/2*dt^2 0
               0 0 1/2*dt^2
               dt 0 0
               0 dt 0
               0 0 dt];

        S=diag([prec.sysNoise^2;prec.sysNoise^2;prec.sysNoise^2]);

        Q_k=N_k*S*N_k';

        % Praediktionsschritt
        x_k_pred=Phi_k*x_k_hat;         % FIXME x_k_hat{i-1}
        P_k_pred=Phi_k*P_k*Phi_k'+Q_k;  % FIXME P_k{i-1}

        % Berechnung der Kalmangewichtsmatrix
        K_k=P_k_pred*H_k'*(H_k*P_k_pred*H_k'+R_k)^(-1);

        % Beobachtungsupdate
        x_k_hat=x_k_pred+K_k*(z_k-H_k*x_k_pred);
        P_k=(eye(length(K_k))-K_k*H_k)*P_k_pred;

end