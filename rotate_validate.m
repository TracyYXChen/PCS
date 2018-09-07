%% -----Info-----
%rotate matrix
%Yuexi (Tracy) Chen
%August 14, 2018
%% -----hyperparams-----
numbat_file = 'data/tensor_numbat.txt';
pcs_exp_file = 'data/pcs_exp.txt';
pdb_file = 'data/1d3z.pdb';
exp_numbat_file = 'data/pcs_exp_numbat.txt';
pdb_model = 1;
which_chi = 'xx';
%numbat para_center
para_center = [56.611 -93.464 -10.031]*10^-10;
%% -----read numbat-----
fid = fopen(numbat_file);
data=textscan(fid,'%s %f %f','HeaderLines',1,'delimiter',' ');
fclose(fid);
chi_ax = data{2}(1);
chi_rh = data{2}(2);
[alpha,beta,gamma] = deal(data{2}(6),data{2}(7),data{2}(8));
%solve chi_xx, chi_yy, chi_zz
chi_zz = 2/3 * chi_ax;
chi_yy = -(chi_rh + chi_zz)/2;
chi_xx = -(chi_yy + chi_zz);
%% -----create rotation matrix-----
%normal basic rotation matrix
Rx = [1 0 0; 0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
Ry = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
Rz = [cosd(gamma) -sind(gamma) 0;sind(gamma) cosd(gamma) 0;0 0 1];
%Here our rotation matrix is z-y-z
Rx = [cosd(alpha) -sind(alpha) 0;sind(alpha) cosd(alpha) 0;0 0 1];
%Rot_mat = Rz*Ry*Rx;
Rot_mat = Rx*Ry*Rz;
%% -----rotate operation-----
diag_chi = diag([chi_xx, chi_yy,chi_zz]);
%we have diag_mat = C.T * A * C
%A = inv(C.T)*diag_mat*inv(C)
%inv(C.T) = C
prev_chi = Rot_mat* diag_chi * Rot_mat';
[prev_chi_xx, prev_chi_yy,prev_chi_zz] = deal(prev_chi(1,1),prev_chi(2,2),prev_chi(3,3));
[prev_chi_xy, prev_chi_xz,prev_chi_yz] = deal(prev_chi(1,2), prev_chi(1,3),prev_chi(2,3));
fprintf('now chi is')
now_chi = [prev_chi_xx, prev_chi_xy, prev_chi_xz, prev_chi_yy, prev_chi_yz, prev_chi_zz];
if which_chi == 'xx'
    chi_mat = [now_chi(2),now_chi(3),now_chi(4), now_chi(5), now_chi(6)]'.*10^-32;
elseif which_chi == 'yy'
    chi_mat = [now_chi(1),now_chi(2),now_chi(3), now_chi(5), now_chi(6)]'.*10^-32;
elseif which_chi == 'zz'
    chi_mat = [now_chi(1),now_chi(2),now_chi(3), now_chi(4), now_chi(5)]'.*10^-32;
else
    fprintf('only xx,yy,zz are supported')
end
%% -----construct A mat-----
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
%here we don't use pcs_exp but use pcs_calc in the future, but just use the residue number
%to extract coordinates
num_res = length(pcs_exp);
pcs_exp = pcs_exp;
A = zeros(num_res,5);
pdb_x = pdb_coor(:,1);
x = pdb_x - para_center(1);
pdb_y = pdb_coor(:,2);
y = pdb_y - para_center(2);
pdb_z = pdb_coor(:,3);
z = pdb_z - para_center(3);
r_sqr = x.^2 + y.^2 +z.^2;
if which_chi == 'zz'
    A(:,1)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [x.^2 - z.^2];
    A(:,2)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*x.*y];
    A(:,3)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*x.*z];
    A(:,4)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [y.^2 - z.^2];
    A(:,5)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* [2.*y.*z];
elseif which_chi == 'xx'
    A(:,1)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*y); 
    A(:,2)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*x.*z); 
    A(:,3)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (y.^2 - x.^2); 
    A(:,4)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (2*y.*z); 
    A(:,5)=(1./r_sqr.^2.5 * 1/(4 * pi)) .* (z.^2 - x.^2); 
elseif which_chi == 'yy'
    fprintf('underconstruction')
else
    fprintf('which_chi could only be xx,yy or zz, others are not supported');
end        
%% -----calculate PCS-----
cond(A)
fprintf('Here is the pcs calculated\n');
cal_pcs = A * chi_mat * 10^6;
%compare with numbat
tmp = dlmread(exp_numbat_file);
fprintf('Here is the pcs returned by Numbat\n');
numbat_pcs = tmp(:,2);
scatter(cal_pcs,pcs_exp);
xlabel("calc PCS/ppm");
ylabel("exp PCS/ppm");
title("My reproduction");




