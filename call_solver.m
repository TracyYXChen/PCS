%define parameter
pdb_coor = pdb_coor;
pcs_exp = pcs_exp;
pcs_file = 'data/pcs_exp_pred.txt';
%call function
guess = [0,0,0];
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
[position, Chi2]=fminsearch(@(guess) pcs_solver(guess,pdb_coor,pcs_exp,pcs_file),guess,options);
%compare with numbat
exp_numbat = dlmread('data/pcs_exp_numbat.txt');
numbat_chi2 = sum((exp_numbat(:,1) - exp_numbat(:,2)).^2);
fprintf('Our chi square is %f, and Chi square of numbat is %f', Chi2, numbat_chi2)
%arguments:
%guess: 
%fit_ratio: experimental data to fit
%fit_coord: protein atom coordinates
%R2dia:
%Htime: