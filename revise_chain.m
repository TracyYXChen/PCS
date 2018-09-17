%read 1d3z.pdb and convert all chainID as 'A_1d3z', 'B_1d3z'
%read mnodes*.pdb and convert all chainID as 'm_A','m_B'
pdb1d3z = 'data/1d3z.pdb';
A_1d3z= 'data/a_1d.pdb';
B_1d3z= 'data/b_1d.pdb';
mnode = 'data/mnodes001.pdb';
di_mnode = 'data/di_mn.pdb';
which = 'A1';
%'B1','dimer'
%mat2pdb(pdb_data,A_1d3z,which);
%pdb_data = pdb2mat(pdb1d3z,1);
%mat2pdb(pdb_data,B_1d3z,'B1');
%pdb_data = pdb2mat(mnode,1);
%mat2pdb(pdb_data,di_mnode,'dimer');
test = 'data/new_1.pdb';
new_pdb = pdb2mat(test);

num_select = length(new_pdb.X);
num_all = length(new_pdb.resNum);
pdb_coor = zeros(72,3);
%create an array to store coordinates of selected atoms
for ii=2: 73
    for jj=1: num_all
      if (new_pdb.resNum(jj)==ii) && isequal(new_pdb.atomName(jj),{'H'})...
              && isequal(new_pdb.chainID(jj),{'y'})
          pdb_coor(ii-1,1) = new_pdb.X(jj)*10^-10; 
          pdb_coor(ii-1,2) = new_pdb.Y(jj)*10^-10;
          pdb_coor(ii-1,3) = new_pdb.Z(jj)*10^-10; 
      end
    end
end
pdb_coor
