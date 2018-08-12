%----------------------------------
%compare numbat tensor with ours
%Yuexi (Tracy) Chen
%August 12, 2018

function [diff_ax,diff_rh, ax_error,rh_error] = read_tensor(numbat_file, chi_file,which_chi)
%read numbat_file
fid = fopen(numbat_file);
data=textscan(fid,'%s %f %f','HeaderLines',1,'delimiter',' ');
fclose(fid);
chi_ax = data{2}(1);
chi_rh = data{2}(2);
[x,y,z] = deal(data{2}(3),data{2}(4),data{2}(5));
[alpha,beta,gamma] = deal(data{2}(6),data{2}(7),data{2}(8));
fid = fopen(chi_file);
my_chi=textscan(fid,'%f %f %f %f %f','delimiter',',');
fclose(fid);
%calculate our chi_ax and chi_rh
%numbat use 10^-32m^-3, here we already have 10^-30m^3, so they need to be
%divided by 100
%ax = zz - (xx+yy)/2
%rh = xx - yy
if which_chi == 'zz'
    my_chiax = (-my_chi{1}-my_chi{4}-(my_chi{1}+my_chi{4})/2)/100;
    my_chirh = (my_chi{1}-my_chi{4})/100;
elseif which_chi == 'xx'
    my_chiax = (my_chi{5}-(-my_chi{3} - my_chi{5} + my_chi{3})/2)/100;
    my_chirh = (-my_chi{3}-my_chi{5} - my_chi{3})/100;
elseif which_chi == 'yy'
    my_chiax = (my_chi{5}-(my_chi{1} - my_chi{5} - my_chi{1})/2)/100;
    my_chirh = (my_chi{1}-(-my_chi{5}-my_chi{1}))/100;
else
    fprintf('only xx, yy, zz are supported types of which_chi')
end
diff_ax = my_chiax - chi_ax;
diff_rh = my_chirh - chi_rh;
ax_error = data{3}(1);
rh_error = data{3}(2);
end

