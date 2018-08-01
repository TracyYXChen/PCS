%call function
guess = [0,0,0];
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
%[position, Chi2]=fminsearch('pcs_solver',guess,options,fit_ratio,fit_coord,factor,R2dia,Htime);
[position, Chi2]=fminsearch('pcs_solver',guess,options);
%arguments:
%guess: 
%fit_ratio: experimental data to fit
%fit_coord: protein atom coordinates
%R2dia:
%Htime: