function [ y ] = func_detectDataGaps( x,t,max_dt )
%FUNC_DETECTDATAGAPS Summary of this function goes here
%   Detailed explanation goes here

    diff_vec=[0 diff(t)'];
    
    j=1;
    
    for i=1:length(diff_vec)
       if diff_vec(i) > max_dt
           temp(1:3,j:j+diff_vec(i)-2) = ...
              [NaN(3,diff_vec(i)-1)];
           j=j+diff_vec(i)-1;
       else
           temp(:,j) = x(:,i);
           j=j+1;
       end
    end
    
    y=temp;
    
end

