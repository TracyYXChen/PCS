%----------------------------------
%Solve Pseudocontact Shift(PCS) tensors 
%Yuexi (Tracy) Chen
%July 31, 2018

function Chi_2 = pcs_solver(guess, pdb_coor, pcs_exp, pcs_file)
%guess: initialization of paramagnetic center
%pbd_coor: pdb coordinates
%pcs_exp: pcs experiment
%pcs_file: file to store exp and predicted pcs for comparison
%construct A matrix
para_center = guess;
num_res = length(pcs_exp);
A = zeros(num_res,5);
x = pdb_coor(:,1);
y = pdb_coor(:,2);
z = pdb_coor(:,3);
for ii = 1:num_res
    r_sqr = (para_center(1) - x(ii))^2 + (para_center(2) - y(ii))^2 +(para_center(3) - z(ii))^2;
    A(ii,:) =(1/r_sqr^2.5 * 1/4 * pi) .* [x(ii)^2 - z(ii)^2, 2*x(ii)*y(ii), 2*x(ii)*z(ii), y(ii)^2 - z(ii)^2, 2*y(ii)*z(ii)];
end
%SVD 
[U, S, V] = svd(A);
%x = A-1 * PCS_exp
%S is not square, use pinv
A_inv = V * pinv(S) * U'; 
x =  A_inv * pcs_exp;
pcs_calc = A * x;
dlmwrite(pcs_file, [pcs_exp, pcs_calc]);
Chi_2 = sum((pcs_calc - pcs_exp).^2);
return


