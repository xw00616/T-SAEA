function [x_pca, M]=PCAMyself(x)

[COEFF, SCORE, latent]=pca(x);
% answer=0;s1=0;s2=sum(latent);count=1;
% while answer < 0.95
%     s1=s1+latent(count);
%     answer=s1/s2;
%     count=count+1;
% end
M=COEFF(:,1:0.5*end);
x_pca=x*M;