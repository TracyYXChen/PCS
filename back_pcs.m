%----------------------------------
%Recover Pseudocontact Shift(PCS) tensors 
%Yuexi (Tracy) Chen
%Aug 13, 2018
%----------------------------------
%--------hyperparameters-----------
which_chi = 'xx';
which_method = 'diag';
chi_file = 'data/chi_file.txt';
pcs_exp_file = 'data/pcs_exp.txt';
pcs_exp_pred_file = 'data/pcs_exp_pred.txt';
pdb_file = 'data/1d3z.pdb';
after_pcs_file = 'data/after_pcs.txt';
numbat_file = 'data/tensor_numbat.txt';
pdb_model = 1;
para_center = [70.8092,-85.3409,-14.6403];
%--------build chi_matrix----------
fid=fopen(chi_file);
data=textscan(fid,'%f %f %f %f %f %f','delimiter',',');
fclose(fid);
if which_chi == 'xx'
    chi_mat = [data{2},data{3},data{4},data{5},data{6}]';
elseif which_chi == 'yy'
    chi_mat = [data{1},data{2},data{3},data{5},data{6}]';
elseif which_chi == 'zz'
    chi_mat = [data{1},data{2},data{3},data{4},data{5}]';
else
    fprintf('only xx,yy,zz are supported')
end
%--------build A matrix----------
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
num_res = length(pcs_exp);
A = zeros(num_res,5);
x = pdb_coor(:,1);
y = pdb_coor(:,2);
z = pdb_coor(:,3);
for ii = 1:num_res
    r_sqr = (para_center(1) - x(ii))^2 + (para_center(2) - y(ii))^2 +(para_center(3) - z(ii))^2;
    if which_chi == 'zz'
        A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [x(ii)^2 - z(ii)^2, 2*x(ii)*y(ii), 2*x(ii)*z(ii), y(ii)^2 - z(ii)^2, 2*y(ii)*z(ii)];
    elseif which_chi == 'xx'
        A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [2*x(ii)*y(ii), 2*x(ii)*z(ii), y(ii)^2 - x(ii)^2, 2*y(ii)*z(ii), z(ii)^2 - x(ii)^2]; 
    elseif which_chi == 'yy'
         A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [x(ii)^2 - y(ii)^2, 2*x(ii)*y(ii),2*x(ii)*z(ii),2*y(ii)*z(ii), z(ii)^2 - y(ii)^2];
    else
        fprintf('which_chi could only be xx,yy or zz, others are not supported');
    end
end
%--------calculate PCS----------
pcs_calc = A * chi_mat ./100;
dlmwrite(after_pcs_file, pcs_calc);
tmp = dlmread(pcs_exp_pred_file);
pcs_prev_calc = tmp(:,2);
fprintf('Are pcs_calc and pcs_rev the same\n');
if pcs_calc == pcs_prev_calc
    fprinf('Yes\n');
else fprintf('No\n');
end
%---------compare position and delta_chi------
guess = [56,-93,-10];
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
fprintf('Start searching process...\n')
[position, Chi2]=fminsearch(@(guess) pcs_solver(guess,pdb_coor,pcs_calc,after_pcs_file, chi_file, which_chi),guess,options);
fprintf('Search finished.\n')
fprintf('Now the position is\n')
position
%compare with numbat tensors
fprintf('Chi-square eliminating type is %s \n',which_chi);
read_tensor(numbat_file, chi_file, which_method);


