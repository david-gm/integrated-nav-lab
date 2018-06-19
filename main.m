%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           lab            %
%   Inertial Navigation    %
%                          %
%        Eva BAUER         %
%        #0830154          %
%       David GMEINDL      %
%        #0831307          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all;

%% load data
pathStaticINS='Daten/INS_static.txt';
pathKinematicINS='Daten/INS_kinematic.txt';
pathGPS='Daten/GPS.txt';
if ~exist('data','var')
    staticData=load(pathStaticINS);
    kinematicData=load(pathKinematicINS);
    gpsData = importfileGPS(pathGPS);
    load('IMUcentreCoords');
    load('constants');
    load('IMUStrapdownConst');
    kalman.in = load('precisions');
end

%% compute initial alignment
[initialAlignment.r,initialAlignment.p,initialAlignment.y,initialAlignment.statVals]=...
    func_InitialAlignment(staticData,IMUcentreCoords,constants);

figure
subplot(3,1,1)
plot(initialAlignment.r,'r');hold on
title('Roll');ylabel('[°]')
subplot(3,1,2)
plot(initialAlignment.p,'g')
title('Pitch');ylabel('[°]')
subplot(3,1,3)
plot(initialAlignment.y,'b')
title('Yaw');ylabel('[°]')

mtit('Time Series of Euler Angles',...
    'fontsize',10,'xoff',0,'yoff',.025,'FontWeight', 'bold');
        


%% compute strapdown algorithm
[imu.x_e,imu.x_e_spherical,imu.t_imu]=...
    func_StrapdownAlgorithm(kinematicData,IMUStrapdownConst,constants);

% difference Vector
% TODO Projections.cartesianToLLF(dX,dY,dZ,phi,lam)
dX = imu.x_e(1,:)' - imu.x_e(1,1);
dY = imu.x_e(2,:)' - imu.x_e(2,1);
dZ = imu.x_e(3,:)' - imu.x_e(3,1);
% ---------------------------------------
% ?????
% ---------------------------------------
%%%imu.x_e_L = Projections.getLocalCoords(imu.x_e, imu.x_e_spherical);
imu.x_e_L = Projections.cartesianToLLF(dX, dY, dZ, imu.x_e_spherical(1,1), imu.x_e_spherical(2,1))';

%% Kalmanfiltering
    for i=1:length(gpsData)
        [a_wgs84, b_wgs84] = Projections.getWGS84Axes();
                
        kalman.in.x_gps(:,i) = Projections.ellipticToCartesianX(a_wgs84, b_wgs84, ...
            [deg2rad(gpsData(i,2:end-1)),gpsData(i,end)])';
        %%%kalman.in.x_gps(:,i) = Projections.getCartesianCoords(...
        %%%[deg2rad(gpsData(i,2:end-1)),gpsData(i,end)],constants.a,constants.f);
    end

    kalman.init.P_initial=eye(6);
   
    
%     kalman.in.x_imu = [[0 0 0];diff(imu.x_e(:,25:50:end)')];
    kalman.in.x_imu = [imu.x_e(:,1250:end)]; % first 25 sec omitted in KF
    kalman.init.x_initial = [   kalman.in.x_gps(:,1)
                                0;0;0];
     
    kalman.in.dt=1;
    
    % find data gaps
    kalman.in.x_gps=func_detectDataGaps(kalman.in.x_gps,gpsData(:,1),kalman.in.dt);
    kalman.in.x_imu=func_detectDataGaps(kalman.in.x_imu,imu.t_imu(1249:end)',kalman.in.dt);
    kalman.plot.gpsDataNoGapsSpherical = func_detectDataGaps(...
        [deg2rad(gpsData(:,2:3))';gpsData(:,4)'],gpsData(:,1),kalman.in.dt);

    kalman.in.x_imu=kalman.in.x_imu';
    
    kalman.in.dt=1/50;
    % kalman filter
    [ kalman.out.x_k_hat ] = func_kalmanFiltering( ...
        kalman.init.x_initial, ...
        kalman.init.P_initial, ...
        kalman.in.x_gps', ...
        kalman.in.x_imu(1:end,:), ...   %FIXME: x_imu(2:end,:)?? should be all or not??
        kalman.in.precisions, ...
        kalman.in.dt );
    
    kalman.out.x_k_hat = [kalman.init.x_initial,kalman.out.x_k_hat];
    
    for i=1:length(kalman.in.x_imu)
        
        [a_wgs84, b_wgs84] = Projections.getWGS84Axes();
        kalman.out.x_k_hat_spherical(:,i) = Projections.cartesianToEllipticX(...
            a_wgs84, b_wgs84, kalman.out.x_k_hat{i}(1:3)', false);
        %%%kalman.out.x_k_hat_spherical(:,i) = Projections.getSphericalCoords(...
        %%%    kalman.out.x_k_hat{i}(1:3),constants.a,constants.f);
    end
    
    % get local coords
    temp = cell2mat(kalman.out.x_k_hat);
    dX = temp(1,:)' - temp(1,1);
    dY = temp(2,:)' - temp(2,1);
    dZ = temp(3,:)' - temp(3,1);
    kalman.out.x_k_hat_L = Projections.cartesianToLLF(dX, dY, dZ, ...
        kalman.out.x_k_hat_spherical(1,1), kalman.out.x_k_hat_spherical(2,1))';

    %%%kalman.out.x_k_hat_L = Projections.getLocalCoords(cell2mat(kalman.out.x_k_hat), ...
    %%%    kalman.out.x_k_hat_spherical);    
            
    dX = kalman.in.x_gps(1,:)' - kalman.in.x_gps(1,1);
    dY = kalman.in.x_gps(2,:)' - kalman.in.x_gps(2,1);
    dZ = kalman.in.x_gps(3,:)' - kalman.in.x_gps(3,1);
    kalman.plot.x_gps_L = Projections.cartesianToLLF(dX, dY, dZ, ...
        kalman.plot.gpsDataNoGapsSpherical(1,1), kalman.plot.gpsDataNoGapsSpherical(2,1))';
    
    %%%kalman.plot.x_gps_L = Projections.getLocalCoords(kalman.in.x_gps, ...
    %%%    kalman.plot.gpsDataNoGapsSpherical);   
        
    % plots
    titles = {'Spherical','Local Level','Heights','Kalman Filter'};
    % FIXME: IMU spherical offset
    func_plotKalman( gpsData, imu, kalman, titles );
    
%% loosly coupled kalman filter    

    % compute strapdown algorithm
            % additional inputs:
            %   1. sampleTimeImuForKalman = 50
            %   2. kalman.in.x_gps       --> should include NaN @ data gaps
            %   3. kalman.init.P_initial
            %   4. kalman.in.precisions
            %   5. kalman.in.dt
            sampleTimeImuForKalman = 50;
            looslyCoupledKF.in.x_gps = kalman.in.x_gps;
            looslyCoupledKF.init.P_initial = kalman.init.P_initial;
            kalman.in.precisions.gps=0.002;
            looslyCoupledKF.in.precisions = kalman.in.precisions;
            looslyCoupledKF.in.dt = kalman.in.dt/50;
            

    [looslyCoupledKF.out.x_e, looslyCoupledKF.out.x_e_spherical, looslyCoupledKF.t_imu]=...
        func_StrapdownAlgorithmCoupled(...
            kinematicData,...
            IMUStrapdownConst,...
            constants,...
            looslyCoupledKF,...
            sampleTimeImuForKalman);


    looslyCoupledKF.plot.gpsDataNoGapsSpherical = kalman.plot.gpsDataNoGapsSpherical;

    % get local coords
    dX = looslyCoupledKF.out.x_e(1,:)' - looslyCoupledKF.out.x_e(1,1);
    dY = looslyCoupledKF.out.x_e(2,:)' - looslyCoupledKF.out.x_e(2,1);
    dZ = looslyCoupledKF.out.x_e(3,:)' - looslyCoupledKF.out.x_e(3,1);
    looslyCoupledKF.out.x_e_L = Projections.cartesianToLLF(dX, dY, dZ, ...
        looslyCoupledKF.out.x_e_spherical(1,1), looslyCoupledKF.out.x_e_spherical(2,1))';
    %%%looslyCoupledKF.out.x_e_L = Projections.getLocalCoords(looslyCoupledKF.out.x_e, ...
    %%%    looslyCoupledKF.out.x_e_spherical);
    dX = looslyCoupledKF.in.x_gps(1,:)' - looslyCoupledKF.in.x_gps(1,1);
    dY = looslyCoupledKF.in.x_gps(2,:)' - looslyCoupledKF.in.x_gps(2,1);
    dZ = looslyCoupledKF.in.x_gps(3,:)' - looslyCoupledKF.in.x_gps(3,1);
    looslyCoupledKF.plot.x_gps_L = Projections.cartesianToLLF(dX, dY, dZ, ...
        looslyCoupledKF.plot.gpsDataNoGapsSpherical(1,1), looslyCoupledKF.plot.gpsDataNoGapsSpherical(2,1))';
    %%%looslyCoupledKF.plot.x_gps_L = Projections.getLocalCoords(looslyCoupledKF.in.x_gps, ...
    %%%    looslyCoupledKF.plot.gpsDataNoGapsSpherical);   

    % plots
    titles = {'Spherical','Local Level','Heights','Loosly Coupled Kalman Filter'};
    % FIXME: IMU spherical offset
    func_plotKalman( gpsData, imu, 0, titles, looslyCoupledKF );    
    
    %savefigs()