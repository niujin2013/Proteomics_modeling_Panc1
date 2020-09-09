%% ciap 
global time y time_s ypred ; 
time = [6 24 48 72 ] ; % repmat(time',1,3); 
% cIAP1 
options = statset('Display' , 'iter' );
ciap = xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cIAP\cIAP.xlsx', 'B3:B14'); 
y = ciap ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for cIAP 
% format short
beta0 = [1 0.941 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitciap',...
beta0 , options);
% visual inspection 

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitciap',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits',out)
global beta_prev ; 
beta_prev  = beta;   
fig = plotfig('cIAP1') 

% save everything 
saveas(fig,'fitciap')
saveas(fig,'fitciap.png')

CV = (diag( sqrt(diag(CovB)).*100 ./beta ) )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;

show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'ciap '); 

%% pNFKB/NFKB
pp65=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\p65nfkb\pnfkbnfkb.xlsx', 'B3:B14'); 
y =pp65; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [ -0.5 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitpp65',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitpp65',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits pp65',out)
fig = plotfig('pp65')  
% save everything 
saveas(fig,'fitpp65')
saveas(fig,'fitpp65.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'pp65 '); 
beta_prev  = [beta_prev beta];  

%% ciap  +  pp65 
y = [ciap; pp65]; 
yp= reshape(y, 4,(size(y,1)/4) ) ;  % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [  0.6084    0.9560  -0.2799 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitciappp65',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitciappp65',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [repmat(time_s,1,(size(y,1)/4)) ; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred ciap&pp65',out)
proteins = ["cIAP1", "pNFkB/NFkB"] ; 
k=0; 


CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
beta_prev  = beta; 
xlswrite('parameter seq',show, 'ciap_pp65 '); 
%% BAX <-pp65 
bax =  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\BAX\BAX.xlsx', 'B3:B14'); 
y =bax ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% fit non linear regression  
beta0 = [0.1 0.5 0.5 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitbax',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitbax',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]
xlswrite('Pred w limits pp65 bax',out)

fig= plotfig('BAX') 

% save everything 
saveas(fig,'fitBAX')
saveas(fig,'fitBAX.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
beta_prev  = [beta_prev beta];  
xlswrite('parameter seq',show, 'bax '); 

%% IRAK4
IRAK4=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\IRAK4\IRAK4.xlsx', 'B3:B14'); 
y = IRAK4 ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [ 0.1  0.2 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitirak4',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitirak4',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits IRAK4',out)
fig= plotfig('IRAK4')  
% save everything 
saveas(fig,'fitIRAK4')
saveas(fig,'fitIRAK4.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'IRAK4'); 
beta_prev  = [beta_prev beta];
%% pJNK 
pJNK=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\pJNK\pJNKJNK.xlsx', 'B3:B14'); 
y = pJNK ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [ 0.1  ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitpjnk',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitpjnk',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits pjnk',out)

fig = plotfig('pJNK'); 
% save everything 
saveas(fig,'fitppjnk')
saveas(fig,'fitpjnk.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'pjnk'); 
beta_prev  = [beta_prev beta];  
%% VDAC1 <-pJNK
beta_prev =0.008263; 
VDAC1=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\VDAC1\VDAC1.xlsx', 'B3:B14'); 
y = VDAC1 ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for VDAC1
beta0 = [ 1 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitVDAC1',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitVDAC1',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits VDAC1',out)

fig = plotfig('VDAC1'); 
% save everything 
saveas(fig,'fitpVDAC1')
saveas(fig,'fitVDAC1.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'VDAC1'); 
beta_prev  = [beta_prev beta];  

%% pSTAT3 <-IRAK4  |- pJNK 
beta_prev  = [0.1270 0.96 0.008263 ]; 
pSTAT3=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\pSTAT3\pSTAT3.xlsx', 'B3:B14'); 
y = pSTAT3 ; 
yp= reshape(y, 4,3) % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [ 0.1 0.2]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitpSTAT3',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitpSTAT3',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits pSTAT3',out)
fig = plotfig('pSTAT3'); 

% save everything 
saveas(fig,'fitpSTAT3')
saveas(fig,'fitpSTAT3.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'pSTAT3'); 


%% BCL2 <-pp65&pSTAT3 & |-VDAC
beta_prev  = forBCL2 ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input parameters from xls file; 
BCL2=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\BCL2\BCL2.xlsx', 'B3:B14'); 
y = BCL2 ; 
yp= reshape(y, 4,3) ;  % for plot purpose
% nonlinear regression for pNFKB/NFKB 
beta0 = [  1.1 1.1 -3 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitBCL2',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitBCL2',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits BCL2',out)
fig  =plotfig('Bcl2'); 

% save everything 
saveas(fig,'fitBCL2')
saveas(fig,'fitBCL2.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'BCL2'); 

%% ELYS <- PTX 
ELYS=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\ELYS\ELYS.xlsx', 'B3:B14'); 
y = ELYS ; 
yp= reshape(y, 4,3) ;  % for plot purpose

% nonlinear regression for pNFKB/NFKB 
beta0 = [  0.1 0.05 1 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitELYS',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitELYS',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits ELYS',out)
fig =plotfig('ELYS'); 

saveas(fig,'fitELYS')
saveas(fig,'fitELYS.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'ELYS'); 

%% ASPP2 <- ELYS
beta_prev = beta; 
ASPP2=  xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\ASPP2\ASPP2.xlsx', 'B3:B14'); 
y = ASPP2 ; 
yp= reshape(y, 4,3) ;  % for plot purpose

% nonlinear regression for pNFKB/NFKB 
beta0 = [ 1 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, y', 'ynewfitASPP2',beta0 ,options);

time_s = 0:0.1:72; 
[ypred,delta] = nlpredci('ynewfitASPP2',time_s , beta,R,'Covar',CovB,'MSE',MSE); 
lower = ypred - delta;
upper = ypred + delta;
out = [[ time_s time_s time_s  ]; lower; ypred; upper ]'; 
Names = ["time_s" "Lower" "Pred" "Upper"  ]; 
out = [Names ; out ]; 
xlswrite('Pred w limits ASPP2',out)
fig = plotfig('ASPP2') ; 

saveas(fig,'fitASPP2')
saveas(fig,'fitASPP2.png')

CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter seq',show, 'ASPP2');



%% fit with cycle
beta_prev = forcasp3 ; 
time = xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cycle\simB-P-BP.xlsx', 'A2:A999'); 
cycle_B = xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cycle\simB-P-BP.xlsx', 'B2:G999'); 
cycle_P = xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cycle\simB-P-BP.xlsx', 'I2:N999'); 
cycle_BP= xlsread('C:\Users\jinniu\Documents\MATLAB\new panc1 proteins\cycle\simB-P-BP.xlsx', 'P2:U999');
cell = [change(cycle_B) change(cycle_P) change(cycle_BP) ]; 
yp= cycle ; 
cell_y = reshape(cell  ,size(cell ,1)*size(cell ,2) , 1) ;

% nonlinear regression:  ELYS_MA BAX_AP  kapm0 cIAP_AP  ktau 
beta0 = [0.02 -0.05 0.12 -0.5 0.05 ]; 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(time, cell_y , 'ynewfitcycle',beta0 ,...
    options  );
param = beta;
CV = (diag( sqrt(diag(CovB)).*100 ./abs(beta))  )';
pars_ci = nlparci(beta,R,'jacobian',J,'alpha',0.05 ) ;
show = [ beta'  CV' pars_ci  ]
xlswrite('parameter cycle',show, 'cycle');

ynewplotcycle

%% simulate/extrapolate other doses in the cell viability assay 
external = timecourseviability ; 


