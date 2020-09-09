% with cycle 
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
function y = ynewfitcycle(param, time  ) 
  % Define beta and delta as globals so we can pass the values into the 
  % ‘IDR1’ function.
 global beta_prev IC_1 IC_2 IC_3 IC_4 IC_5 IC_6 ; 
 beta_prev = xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cIAP\parameter seq.xls', 'B2:B20');
 global ELYS_MA ASPP2_AP  kapm0 cIAP_AP ktau; 
 
 ELYS_MA = param(1) ; 
 ASPP2_AP  = param(2) ;  
 kapm0 = param(3) ; 
 cIAP_AP = param(4) ; 
 ktau = param(5) ;  
 
        IC_1    =   0.2617; 
        IC_2    =   0.1468; 
        IC_3    =   0.489; 
        IC_4 = 0; 
        IC_5 = 0; 
        IC_6 = 0.0903;  
  ICC = [IC_1 IC_2 IC_3 IC_4 IC_5 IC_6]    ; 
  i = 40; 
  ICS = [ones(1,i) ICC ICC ICC]; 
  % Solve our model.  
  [t, y] = ode15s(@IDRnewcycle, [0 time(end)], ICS) ; %ODE(dF, [T], [ICS])
   
 

forB = [y(:, 41)  y(:, 42) y(:, 43)+y(:, 44) y(:, 45) y(:, 46) ]; 
forP = [y(:, 47)  y(:, 48) y(:, 49)+y(:, 50) y(:, 51) y(:, 52) ]; 
forBP= [y(:, 53)  y(:, 54) y(:, 55)+y(:, 56) y(:, 57) y(:, 58) ]; 
forall = [forB forP forBP ]; 
forallint =  interp1(t, forall, time ) ;
y = reshape(forallint, size(forallint,1 )*size(forallint,2 ) ,1)  ; 


          
