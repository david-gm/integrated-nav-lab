function func_plotKalman( gpsData, imu, kalman, titles, varargin )
%FUNC_PLOTKALMAN Summary of this function goes here
%   Detailed explanation goes here

    if nargin == 4  % plot normal kalman filter

        figure
        plot(kalman.plot.x_gps_L(2,:),kalman.plot.x_gps_L(1,:),'x');hold on
        plot(imu.x_e_L(2,:),imu.x_e_L(1,:),'-g')
        plot(kalman.out.x_k_hat_L(2,:),kalman.out.x_k_hat_L(1,:),'-r')
        legend('GPS','IMU','IMU + GPS','location','best')
        xlabel('y [m]');ylabel('x [m]');
        title(titles{2})
        axis equal

        mtit(titles{4},...
             'fontsize',10,'xoff',0,'yoff',.025,'FontWeight', 'bold');
        
        figure
        plot(kalman.plot.gpsDataNoGapsSpherical(3,25:end),'b');hold on;
        plot(imu.x_e_spherical(3,1250:50:end),'g')
        plot(kalman.out.x_k_hat_spherical(3,25:50:end),'r')
        legend('GPS','IMU','IMU+GPS','location','best')
        xlabel('Number of Observations');ylabel('Height [m]');
        title(titles{3})
        
        mtit(titles{4},...
             'fontsize',10,'xoff',0,'yoff',.025,'FontWeight', 'bold');
        
    else    % plot loosly coupled kalman filter
        coupledKalman = varargin{1};
        
        figure
%         subplot(2,1,1)
%         plot(gpsData(:,3),gpsData(:,2),'xb');hold on
%         plot(rad2deg(imu.x_e_spherical(2,:)),rad2deg(imu.x_e_spherical(1,:)),'-g')
%         plot(rad2deg(coupledKalman.out.x_e_spherical(2,:)),rad2deg(coupledKalman.out.x_e_spherical(1,:)),'-r')
%         legend('GPS','IMU','IMU + GPS')
%         title(titles{1})
%         axis equal
% 
%         subplot(2,1,2)
        plot(coupledKalman.plot.x_gps_L(2,:),coupledKalman.plot.x_gps_L(1,:),'x');hold on
        plot(imu.x_e_L(2,:),imu.x_e_L(1,:),'-g')
        plot(coupledKalman.out.x_e_L(2,:),coupledKalman.out.x_e_L(1,:),'-r')
        legend('GPS','IMU','IMU + GPS','location','best')
        xlabel('y [m]');ylabel('x [m]');
        title(titles{2})
        axis equal

        mtit(titles{4},...
             'fontsize',10,'xoff',0,'yoff',.025,'FontWeight', 'bold');
        
        figure
        plot(coupledKalman.plot.gpsDataNoGapsSpherical(3,:),'b');hold on;
        plot(imu.x_e_spherical(3,25:50:end),'g')
        plot(coupledKalman.out.x_e_spherical(3,25:50:end),'r')
        legend('GPS','IMU','IMU+GPS','location','best')
        title(titles{3})
        
        mtit(titles{4},...
             'fontsize',10,'xoff',0,'yoff',.025,'FontWeight', 'bold');
        
    end

end

