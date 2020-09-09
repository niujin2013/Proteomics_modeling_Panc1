% ynew fit.m 
%
% Returns expected infection values given a set of time points and 
%   coefficients 
% Inputs: 
% a – List of free parameter values 
% x – List of time points of interest 
% Output: 
% y – A list of expected value for each time 
%     point in ‘x’. 
function y = ynewfitELYS(param, time  ) 
  % Define beta and delta as globals so we can pass the values into the 
  % ‘IDR1’ function. 
 global kdeg_ELYS ELYS_P fb_selfELYS; % = param() ;   
  kdeg_ELYS = param(1) ;
  ELYS_P = param(2) ; 
  fb_selfELYS= param(3) ; 
  
  i =6 ;   
  % Solve our model.  
  [t, y] = ode15s(@IDRnewELYS , [0 time(end)], ones(1,i  )  ) ; %ODE(dF, [T], [ICS])

  % We will use the built-in interpolation function interp1 to get 
  % values of infected people at that correspond to the proper time 
  % points. 
  
  
  y = log2([ interp1(t, y(:,1) , time ) , ...
             interp1(t, y(:,2) , time ) , ...
             interp1(t, y(:,3) , time ) , ...           
             ]) ;
         
         

%              , ... ,...
%              interp1(t, y(:,16) , x),...
%              interp1(t, y(:,17) , x), ...
%              interp1(t, y(:,18) , x),...
%              interp1(t, y(:,19) , x)  
%              interp1(t, y(:,13) , x),...
%              interp1(t, y(:,14) , x),...
%              interp1(t, y(:,15) , x) ,...
%              interp1(t, y(:,16) , x),...
%              interp1(t, y(:,17) , x), ...
%              interp1(t, y(:,18) , x),...
%              interp1(t, y(:,19) , x)     
