function [r,p,y,statisticalValues] = func_InitialAlignment( data, IMUcentre,constants )
%FUNC_INITIALALIGNMENT Summary of this function goes here
%   Detailed explanation goes here


    t=data(:,1);
    
    %rotation (initial-body frame in body frame)
    omega_ib_b=[data(:,3),data(:,2),-data(:,4)];
    
    %specific force(body frame)
    f_b=[data(:,6),data(:,5),-data(:,7)];
    
    [t,omega_ib_b,f_b]=func_Averaging(t,omega_ib_b,f_b);
    
    % is des richtig? unabhängig von t, daher nur für die initial epoche??
    % d.h. dim(omega)=3x1 und dim(g)=3x1??
    % eigentlich müsste t zumindest in w_e drinnen sein, oda?
    [omega_ie_l,g_l]=func_getModel_omega_g([deg2rad(IMUcentre.sphericalCoords(1:2)),IMUcentre.sphericalCoords(3)],constants);
    
    c_b = cross(-f_b,omega_ib_b);
    c_l = cross(g_l,omega_ie_l);
    F_l_inv = inv([g_l,omega_ie_l,c_l]);
    
    for i=1:length(c_b)
        F_b=[-f_b(i,:)',omega_ib_b(i,:)',c_b(i,:)'];
        R_bl(:,:,i)=F_b*F_l_inv;  
    end
    
    r=rad2deg(squeeze(atan2(R_bl(2,3,:),R_bl(3,3,:))));
    p=rad2deg(squeeze(asin(-R_bl(1,3,:))));
    y=rad2deg(squeeze(atan2(R_bl(1,2,:),R_bl(1,1,:))));

    statisticalValues.mean=mean([r p y]);
    statisticalValues.std=std([r p y]);
end

function [t,omega,f] = func_Averaging(t_,omega_,f_)

    dt=diff(t_);

    sum=0;
    k=1;
    beginningIndex=1;
    
    for i=1:length(dt)
        
        sum=sum+dt(i);
        if sum==1 
            sum=0;
            t(k)=(t_(i+1)+t_(beginningIndex))/2;
            omega(k,:)=mean(omega_(beginningIndex:i+1,:));
            f(k,:)=mean(f_(beginningIndex:i+1,:));
            k=k+1;
            beginningIndex=i+1;
        end
    end
        

end