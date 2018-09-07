%% -----Info-----
%perform experiment
%-------
%outline
%-------
%1. read PDB file and extract coordinates according to residues in PCS file
%2. Initialize coordinates of paramagnetic center
%3. solve PCS by SVD 
%4. calculate Chi^2 between PCS_exp and PCS_theory
%5. reinitialize coor of paramagnetic center
%6. After several iterations, return calculated PCS tensors and coor of paramagnetic center
%% -----hyperparams-----
pcs_exp_file = 'data/pcs_exp.txt';
pdb_file = 'data/1d3z.pdb';
pdb_model = 1;
%1st column is exp, 2nd column is numbat prediction
exp_numbat_file = 'data/pcs_exp_numbat.txt';
%output file: 1st column exp, 2nd column: our prediction value
pcs_exp_pred_file = 'data/pcs_exp_pred.txt';
numbat_file = 'data/tensor_numbat.txt';
our_tensor_file = 'data/tensor_ours.txt';
chi_file = 'data/chi_file.txt';
which_chi = 'xx';
%which_method could only be 'diag' or 'formula'
which_method = 'diag';
%extract coordinates and exp value
%% -----simplex-----
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
%call function
guess = [57 -94 -9]*10^-10;
numbat_posi = [56.611 -93.464 -10.031]*10^-10;
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
fprintf('Start searching process...\n')
[position, Chi2]=fminsearch(@(guess) svd_solver(guess,pdb_coor,pcs_exp,pcs_exp_pred_file, chi_file, which_chi),guess,options);
fprintf('Search finished.\n')
fprintf('Now the position is')
position
fprintf('The distance between prediction and numbat is')
ditance = sqrt(sum((numbat_posi - position).^2))
%compare with numbat tensors
fprintf('Chi-square eliminating type is %s \n',which_chi);
read_tensor(numbat_file, chi_file, which_method);
%compare with numbat pcs results
%exp_numbat = dlmread(exp_numbat_file);
%numbat_chi2 = sum((exp_numbat(:,1) - exp_numbat(:,2)).^2);
%fprintf('Our chi square is %f, and Chi square of numbat is %f', Chi2, numbat_chi2)
