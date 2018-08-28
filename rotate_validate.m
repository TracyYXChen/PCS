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
which_chi = 'zz';
%numbat para_center
para_center = [56.611 -93.464 -10.031];
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
%Rot_mat = Rx*Ry*Rz;
%% -----rotate operation-----
diag_chi = diag([chi_xx, chi_yy,chi_zz]);
%we have diag_mat = C.T * A * C
%A = inv(C.T)*diag_mat*inv(C)
%inv(C.T) = C
prev_chi = Rot_mat * diag_chi * inv(Rot_mat);
[prev_chi_xx, prev_chi_yy,prev_chi_zz] = deal(prev_chi(1,1),prev_chi(2,2),prev_chi(3,3));
[prev_chi_xy, prev_chi_xz,prev_chi_yz] = deal(prev_chi(1,2), prev_chi(1,3),prev_chi(2,3));
fprintf('now chi is')
now_chi = [prev_chi_xx, prev_chi_xy, prev_chi_xz, prev_chi_yy, prev_chi_yz, prev_chi_zz];
if which_chi == 'xx'
    chi_mat = [now_chi(2),now_chi(3),now_chi(4), now_chi(5), now_chi(6)]';
elseif which_chi == 'yy'
    chi_mat = [now_chi(1),now_chi(2),now_chi(3), now_chi(5), now_chi(6)]';
elseif which_chi == 'zz'
    chi_mat = [now_chi(1),now_chi(2),now_chi(3), now_chi(4), now_chi(5)]';
else
    fprintf('only xx,yy,zz are supported')
end
%% -----construct A mat-----
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
%here we don't use pcs_exp but use pcs_calc in the future, but just use the residue number
%to extract coordinates
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
%% -----calculate PCS-----
fprintf('Here is the pcs calculated')
cal_pcs = A * chi_mat
%compare with numbat
tmp = dlmread(exp_numbat_file);
fprintf('Here is the pcs returned by Numbat')
numbat_pcs = tmp(:,2);




