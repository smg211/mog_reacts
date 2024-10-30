function react = pca_react(train, test)
% Perform the basic steps of the PCA reactivations analysis outlined in
% Peyrache et al., 2010 J Comp. Neuro
%
% Use as: 
%       react = pca_react(train, test)
%
% Where:
%       train = Neuron x Time matrix of the data from the template period
%       test = Neuron x Time matrix of the data from the match period
%
% Sandon Griffin (2024)

% Initializie reactivations structure to store output
react = [];

% Compute pairwise cell activity correlation matrix for training
C = (1/size(train, 2))*train*train';
C(isnan(C)) = 0;

% Save the correlation matrix
react.corr_mat = C;

% PCA
[V, D_eig] = eig(C);
[~, ~, s2_explain] = pcacov(C);
lambda = D_eig(logical(eye(size(D_eig,1))));

% Determine significance of PCs
q = size(train,2)/size(train,1); % defined by Marcenko-Pastur distribution
s2 = (std(train, 0, 'all'))^2;
lambda_max = s2*((1+sqrt(1/q))^2); % threshold for eigenvalue of signal PCs vs non-signal PCs, as determined by MP distribution

% Re-order PCs by Eigenvalues
[lambda, idx] = sort(lambda, 'descend');
V = V(:, idx);
s2_explain = s2_explain(idx);

% Identity matrix indicating significance of PCs
lambda_sig = lambda > lambda_max;
pc2proj = find(lambda_sig);

% Determine projectors from template
P = {};
for l = 1:size(V, 1) % for each PC
  P{l} = V(:,l)*V(:,l)';
end

% Apply projectors to data
R = nan(length(pc2proj), size(test, 2));
for p = 1:length(pc2proj) % for each pc
  P_zerodiag = P{pc2proj(p)};
  P_zerodiag(find(eye(size(P_zerodiag, 1)))) = 0;
  R(p, :) = sum(test'*P_zerodiag.*test',2)';
end

% Store output
react.pca.sig_pcs = lambda_sig;
react.pca.eigenvector = V;
react.pca.eigenvalue = lambda;
react.pca.s2_explained = s2_explain;
react.projectors.projectors = P;
react.react = R;
