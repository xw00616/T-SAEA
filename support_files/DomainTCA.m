function X_tar_new=DomainTCA(PopDecS,PopDecT,D)
%% Load data
fts=PopDecS;
fts = fts ./ repmat(sum(fts,2),1,size(fts,2));
Xs = zscore(fts,1);    clear fts
%             Ys = PopObj1;           clear labels
fts= PopDecT;
fts = fts ./ repmat(sum(fts,2),1,size(fts,2));
Xt = zscore(fts,1);     clear fts
%             Yt = PopObj2;            clear labels
%% Set algorithm options
options.gamma = 1.0;
options.lambda = 0.1;
options.kernel_type = 'linear';
options.T = 6;
options.dim = D;
options.mu = 0;
options.mode = 'W-BDA';
%% Run algorithm
[X_src_new,X_tar_new,W] = TCA(Xs,Xt,options);
%%%%%%%
end