%% -----Info-----
%Solve Pseudocontact Shift(PCS) tensors 
%Yuexi (Tracy) Chen
%July 31, 2018
%% -----solve svd -----
function [Chi_2] = svd_solver(guess, pdb_coor, pcs_exp, pcs_file, chi_file, which_chi)
%guess: initialization of paramagnetic center
%pbd_coor: a list of pdb coordinates
%pcs_exp: pcs experiment
%pcs_file: file to store exp and predicted pcs for comparison
%construct A matrix
para_center = guess;
num_res = length(pcs_exp);
A = zeros(num_res,5);
pdb_x = pdb_coor(:,1);
pdb_y = pdb_coor(:,2);
pdb_z = pdb_coor(:,3);
for ii = 1:num_res
    x = pdb_x(ii) - para_center(1);
    y = pdb_y(ii) - para_center(2);
    z = pdb_z(ii) - para_center(3);
    r_sqr = x^2 + y^2 +z^2;
    if which_chi == 'zz'
        A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [x^2 - z^2, 2*x*y, 2*x*z, y^2 - z^2, 2*y*z];
    elseif which_chi == 'xx'
        A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [2*x*y, 2*x*z, y^2 - x^2, 2*y*z, z^2 - x^2]; 
    elseif which_chi == 'yy'
         A(ii,:)=(1/r_sqr^2.5 * 1/(4 * pi)) .* [x^2 - y^2, 2*x*y,2*x*z,2*y*z, z^2 - y^2];
    else
        fprintf('which_chi could only be xx,yy or zz, others are not supported');
    end
end       
%SVD
fprintf('condition number of A is %f\n',cond(A));
%[U, S, V] = svd(A);
%x = A-1 * PCS_exp
%S is not square, use pinv
%A_inv = V * pinv(S) * U';
%x =  A_inv * pcs_exp;
if which_chi == 'zz'
    [chi_xx,chi_xy,chi_xz,chi_yy,chi_yz] = deal(x(1),x(2),(3),x(4),x(5));
    chi_zz = -(chi_xx + chi_yy);
elseif which_chi == 'xx'
    x = 10^-32.*[-1.4025, -3.2879,-9.6487,-2.6452,4.8869];
    %[chi_xy,chi_xz,chi_yy,chi_yz,chi_zz] = deal(x(1),x(2),(3),x(4),x(5));
    %chi_xx = -(chi_yy + chi_zz);
elseif which_chi == 'yy'
    [chi_xx,chi_xy,chi_xz,chi_yz,chi_zz] = deal(x(1),x(2),(3),x(4),x(5));
    chi_yy = -(chi_xx + chi_zz);
else
    fprintf('which_chi could only be xx,yy or zz, others are not supported');
end
pcs_calc = A * x';
%dlmwrite(chi_file, [chi_xx, chi_xy,chi_xz,chi_yy,chi_yz, chi_zz]);
norm_pcs_calc = 10^6 * pcs_calc;
dlmwrite(pcs_file, [norm_pcs_calc]);
Chi_2 = sum((pcs_calc - pcs_exp).^2);
return


