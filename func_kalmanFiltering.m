function [ x_k_hat ] = func_kalmanFiltering( x_initial, P_initial, x_gps, x_imu, prec, dt )
%FUNC_KALMANFILTERING Summary of this function goes here
%   Detailed explanation goes here
    
    x_k_hat{1}= x_initial;
    P_k{1}=P_initial;
%     t=0;
    t=25; % >> eigentlich doch 25 sec???

    i_gps=0;
    
    for i=2:length(x_imu)
        Phi_k=[   1 0 0 dt 0 0
                  0 1 0 0 dt 0
                  0 0 1 0 0 dt
                  0 0 0 1 0 0
                  0 0 0 0 1 0
                  0 0 0 0 0 1];
              
% FIXME: in geschwindigkeit
        sx=(0.5*1.4*10^(-3)*t^2+1/6*3.8*10^(-8)*t^3); 
        
        if mod(i,50)==0
            i_gps=i_gps+1;
        end
        if mod(i,50)==0 && ~isnan(x_gps(i_gps,1))
            H_k=[eye(3),zeros(3);eye(3),zeros(3)];
            R_k=diag([prec.gps^2;prec.gps^2;prec.gps^2;sx^2;sx^2;sx^2]);
            
            % Beobachtungsvektor
            z_k=[ x_gps(i_gps,:)'
                  x_imu(i,:)'];
        else
            H_k=[eye(3),zeros(3)];
            R_k=diag([sx^2;sx^2;sx^2]);
            
            % Beobachtungsvektor
            z_k=x_imu(i,:)';
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
        x_k_pred=Phi_k*x_k_hat{i-1};
        P_k_pred=Phi_k*P_k{i-1}*Phi_k'+Q_k;

        % Berechnung der Kalmangewichtsmatrix
        K_k=P_k_pred*H_k'*(H_k*P_k_pred*H_k'+R_k)^(-1);

        % Beobachtungsupdate
        x_k_hat{i}=x_k_pred+K_k*(z_k-H_k*x_k_pred);
        P_k{i}=(eye(length(K_k))-K_k*H_k)*P_k_pred;
        
        t=t+dt;
    end

end

