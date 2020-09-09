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
function y = ynewfitpSTAT3(param, time  ) 
  % Define beta and delta as globals so we can pass the values into the 
  % ‘IDR1’ function.
 global kdeg_pSTAT3 pSTAT3_pJNK  ;
 
  kdeg_pSTAT3 = param(1) ;
  pSTAT3_pJNK   = param(2) ;

  
  i =9 ;   
  % Solve our model.  
  [t, y] = ode15s(@IDRnewpSTAT3 , [0 time(end)], ones(1,i  )  ) ; %ODE(dF, [T], [ICS])

  % We will use the built-in interpolation function interp1 to get 
  % values of infected people at that correspond to the proper time 
  % points. 
  
  y = log2([ interp1(t, y(:,7) , time ) , ...
             interp1(t, y(:,8) , time ) , ...
             interp1(t, y(:,9) , time ) , ...           
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
