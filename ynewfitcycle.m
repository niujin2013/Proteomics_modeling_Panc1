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
   
%         live_B = y(:,41) +y(:,42) +y(:,43)+ y(:,44)+ y(:,45) ; %  !total live cells :B
%         live_P = y(:,47) +y(:,48) +y(:,49)+ y(:,50)+ y(:,51) ; %  !total live cells :P
%         live_BP= y(:,53) +y(:,54) +y(:,55)+ y(:,56)+ y(:,57) ; %  !total live cells :C
%         
%   cycle_B = y(:, 41:45) ./ live_B .*100 ; 
%   cycle_P = y(:, 47:51) ./ live_P .*100 ;  
%   cycle_BP= y(:, 53:57) ./ live_BP.*100 ;    
%   
%   cycle_B1 = [cycle_B(:, 1:2) cycle_B(:, 3)+ cycle_B(:,4) cycle_B(:, 5)] ; 
%   cycle_P1 = [cycle_P(:, 1:2) cycle_P(:, 3)+ cycle_P(:,4) cycle_P(:, 5)] ;   
%   cycle_BP1= [cycle_BP(:, 1:2) cycle_BP(:, 3)+ cycle_BP(:,4) cycle_BP(:, 5)] ;   
%   
%   apo_B =  y(:, 46) ./ (live_B+ y(:, 46)) *100 ; 
%   apo_P =  y(:, 52) ./ (live_P+ y(:, 52)) *100 ;   
%   apo_BP=  y(:, 58) ./ (live_BP+y(:, 58)) *100 ; 
%   
%   N_B  =(live_B+ y(:, 46)) ; 
%   N_P  =(live_P+ y(:, 52)) ;
%   N_BP =(live_BP+y(:, 58)) ;

%OUT   
%   OUT = [cycle_B1 apo_B N_B cycle_P1 apo_P N_P cycle_BP1 apo_BP N_BP ]; %PERCENTAGES      
%   ALL = interp1(t, y , time ) ;  
%   y =  interp1(t, OUT, time ) ;
%  
%   y =  reshape(y, size(y,1 )*size(y,2 ) ,1)  ; 

forB = [y(:, 41)  y(:, 42) y(:, 43)+y(:, 44) y(:, 45) y(:, 46) ]; 
forP = [y(:, 47)  y(:, 48) y(:, 49)+y(:, 50) y(:, 51) y(:, 52) ]; 
forBP= [y(:, 53)  y(:, 54) y(:, 55)+y(:, 56) y(:, 57) y(:, 58) ]; 
forall = [forB forP forBP ]; 
forallint =  interp1(t, forall, time ) ;
y = reshape(forallint, size(forallint,1 )*size(forallint,2 ) ,1)  ; 


          
