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

%% -----hyper parameters-----
pcs_exp_file = 'data/pcs_exp.txt';
pdb_file = 'data/1d3z.pdb';
pdb_model = 1;
%1st column is exp, 2nd column is numbat prediction
exp_numbat_file = 'data/pcs_exp_numbat.txt';
%output file: 1st column exp, 2nd column: our prediction value
pcs_exp_pred_file = 'data/pcs_exp_pred.txt';
numbat_file = 'data/tensor_numbat.txt';
chi_file = 'data/chi_file.txt';
which_chi = 'xx';
which_method = 'diag';
%% -----preprocess-----
tic
%extract coordinates and exp value
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
%numbat = [56.61	-93.46	-10.03	2.90	-2.08	1.18	7.34	0.5];
position_and_chi = zeros(1,8);
%%only need to initialize 5 vars in chi_tensor
%value ranges from (-2numbat,2numbat)
position_and_chi(1) = 56.61*2 - 4*56.61*rand(1);
position_and_chi(2) = -93.46*2 + 4*93.46*rand(1);
position_and_chi(3) = -10.03*2 + 4*10.03*rand(1);
position_and_chi(4) = 2.9*2 - 4*2.9*rand(1);
position_and_chi(5) = -2.08*2 + 4*2.08*rand(1);
position_and_chi(6) = 1.18*2 - 4*1.18*rand(1);
position_and_chi(7) = 7.34*2 - 4*7.34*rand(1);
position_and_chi(8) = 0.5*4 - 4*0.5*rand(1);
fprintf("initialization is")
position_and_chi
%numbat_posi = [56.611 -93.464 -10.031];
%% -----simplex-----
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'MaxIter',100000);
fprintf('Start searching process...\n')
[position_and_chi, Chi2]=fminsearch(@(guess) guess_solver(guess,pdb_coor,pcs_exp,pcs_exp_pred_file, chi_file, which_chi),guess,options);
fprintf('Search finished.\n')
read_tensor(numbat_file, chi_file,which_method)
fprintf('Now the position is')
position = position_and_chi(1:3)
fprintf('Now the chi_tensor is')
chi_tensor = position_and_chi(4:8)
toc
