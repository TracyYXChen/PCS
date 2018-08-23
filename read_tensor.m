%----------------------------------
%compare numbat tensor with ours
%Yuexi (Tracy) Chen
%August 12, 2018

function read_tensor(numbat_file, chi_file,which_method)
%read numbat_file
fid = fopen(numbat_file);
data=textscan(fid,'%s %f %f','HeaderLines',1,'delimiter',' ');
fclose(fid);
fprintf('Below are ax and rh tensor returned by Numbat\n')
chi_ax = data{2}(1)
chi_rh = data{2}(2)
fprintf('Below are positions returned by Numbat\n')
[x,y,z] = deal(data{2}(3),data{2}(4),data{2}(5));
[alpha,beta,gamma] = deal(data{2}(6),data{2}(7),data{2}(8));
fid = fopen(chi_file);
my_chi=textscan(fid,'%f %f %f %f %f %f','delimiter',',');
fclose(fid); 
[chi_xx,chi_xy,chi_xz,chi_yy,chi_yz,chi_zz] = deal(my_chi{1},my_chi{2},my_chi{3},my_chi{4},my_chi{5},my_chi{6});
fprintf('Here are chi_xx, chi_xy, chi_xz, chi_yy, chi_yz and chi_zz');
[chi_xx,chi_xy,chi_xz,chi_yy,chi_yz,chi_zz]
if strcmp(which_method,'formula')
%Method 1: calculate by formula below
%ax = zz - (xx+yy)/2
%rh = xx - yy
    fprintf('Below are my_ax and my_rh returned by formula')
    my_chiax = (chi_zz-(chi_xx+chi_yy)/2)
    my_chirh = (chi_xx-chi_yy)
%diff_ax = my_chiax - chi_ax;
%diff_rh = my_chirh - chi_rh;
%Method 2: calculate by diag
elseif strcmp(which_method, 'diag')
    my_mat = zeros(3);
    [my_mat(1,1),my_mat(2,2),my_mat(3,3)] = deal(chi_xx,chi_yy,chi_zz);
    [my_mat(1,2),my_mat(2,1)] = deal(chi_xy,chi_xy);
    [my_mat(1,3),my_mat(3,1)] = deal(chi_xz,chi_xz);
    [my_mat(2,3),my_mat(2,3)] = deal(chi_yz,chi_yz);
    fprintf('Below are my_ax and my_rh returned by diag');
    my_eig = eig(my_mat);
    my_chiax = my_eig(3) - (my_eig(1) + my_eig(2))/2
    my_chirh = my_eig(1) - my_eig(2)
else
    fprinf('%s is not supported, only accept diag and formula',which_method)
ax_error = data{3}(1);
rh_error = data{3}(2);
end

