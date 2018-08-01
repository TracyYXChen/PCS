%----------------------------------
%Solve Pseudocontact Shift(PCS) tensors 
%Yuexi Chen
%July 31, 2018
%----------------------------------
%1. read PDB file and extract coordinates + read PCS file
%2. Initialize coordinates of paramagnetic center
%3. solve PCS by SVD 
%4. calculate Chi^2 between PCS_exp and PCS_theory
%5. reinitialize coor of paramagnetic center
%6. After several iterations, return calculated PCS tensors and coor of paramagnetic center

%------------
%read coordinates
%Suppose we have N nucleaus, x, y, z should be (N*1) 
function Chi_2 = pcs_svd(pdb_coor, pcs_exp) 
[x, y, z] = read_coor();
para_center = init();
%construct A matrix
A = zeros(N,5);
for ii = 1:N
    r_sqr = (para_center(1) - x(ii))^2 + (para_center(2) - y(ii))^2 +(para_center(3) - z(ii))^2;
    A(:,ii) =(1/r_sqr^2.5 * 1/4 * pi) .* [x(ii)^2 - z(ii)^2, 2*x(ii)*y(ii), 2*x(ii)*z(ii), y(ii)^2 - z(ii)^2, 2*y(ii)*z(ii)];
end
%SVD 
[U, S, V] = svd(A);
%x = A-1 * PCS_exp
A_inv = V * diag(1./diag(S)) * U'; 
x =  A_inv * pcs_exp;
pcs_calc = A * x;
Chi_2 = sum((pcs_calc - pcs_exp).^2);
return


