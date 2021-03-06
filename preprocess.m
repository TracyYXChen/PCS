%% -----text handler-----
%input: 
%1. pcs_exp, 4 columns: residue number, atom name, pcs_exp, error?
%2. corresponding PDB file   
%output: coordinates of those atoms in the PDB file

%% -----read pcs file -----
function [pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model)
fid = fopen(pcs_exp_file);
data=textscan(fid,'%f %s %f %f','delimiter',' ');
fclose(fid);
res_num = data{1};
atom = data{2};
pcs_exp = data{3};
err = data{4};
%---read pdb file ---
PDBdata = pdb2mat(pdb_file,pdb_model);
fprintf('Finished reading PDB file.\n')
%num of all atoms and selected atoms
num_select = length(res_num);
num_all = length(PDBdata.resNum);
%create an array to store coordinates of selected atoms
pdb_coor = zeros(num_select,3);
for ii=1: num_select
    for jj=1: num_all
      if (res_num(ii) == PDBdata.resNum(jj)) && isequal(atom(ii),PDBdata.atomName(jj))
          pdb_coor(ii,1) = PDBdata.X(jj)*10^-10; 
          pdb_coor(ii,2) = PDBdata.Y(jj)*10^-10;
          pdb_coor(ii,3) = PDBdata.Z(jj)*10^-10; 
      end
    end
end

