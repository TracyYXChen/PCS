%% -----Info-----
%discrete optimization
%Yuexi (Tracy) Chen
%August 23, 2018
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
%% -----put into SVD solver and return smallest-----
[pcs_exp,pdb_coor] = preprocess(pcs_exp_file, pdb_file, pdb_model);
%% -----put into SVD solver-----
best_posi = [0,0,0];
minChi2 = 10000000.0;
minresi = 0;
fprintf('Start iterating all pdb coordinates...')
for ii=1:length(pcs_exp)-1
prev_posi = [pdb_coor(ii,1), pdb_coor(ii,2), pdb_coor(ii,3)];
next_posi = [pdb_coor(ii+1,1), pdb_coor(ii+1,2), pdb_coor(ii+1,3)];
now_posi  = (prev_posi + next_posi)/2;
Chi2 = svd_solver(now_posi,pdb_coor,pcs_exp,pcs_exp_pred_file, chi_file, which_chi)
if Chi2 < minChi2
    minChi2 = Chi2;
    minresi = ii;
    best_posi = now_posi;
end
end
fprintf('iteration finished.')
%% -----results-----
g=sprintf('%f ', best_posi);
fprintf('Now the %dth and %dth residue is the best position\n the coordinates are: %s\n the Chi2 is %f\n',minresi,minresi+1,g,minChi2);
%% -----compare with numbat-----
fprintf('----------------')
numbat_posi = [56.611 -93.464 -10.031];
u=sprintf('%f ', numbat_posi);
tmp_pcs_numbat = dlmread(exp_numbat_file);
numbat_chi2 = sum((tmp_pcs_numbat(:,2) - pcs_exp*10^6).^2);
fprintf('Numbat results\n coordinates of the best position are %s\n the Chi2 is %f\n',u,numbat_chi2);
