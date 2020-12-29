function [Popmean,MSE] = FitPre(PopDec,M,Model)
% Update the offline archive in PICEA-g

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
[N,~]  = size(PopDec);
Popmean = zeros(N,M);
MSE    = zeros(N,M);
% PopmeanCO = zeros(N,M);
% MSECO    = zeros(N,M);
[Popmean(:,1),MSE(:,1)]     = predict(Model{1},PopDec);%   use mean
[Popmean(:,2),MSE(:,2)]     = predict(Model{2},PopDec);
% [PopmeanCO(:,1),MSECO(:,1)] = predict(Model{2},PopDec);
% [PopmeanCO(:,2),MSECO(:,2)] = predict(Model{3},PopDec);
% Popmean(:,2)=mean(PopmeanCO(:,1)+PopmeanCO(:,2),2);
% MSE(:,2)=  mean(MSECO(:,1)+MSECO(:,2),2);
% Popmean(:,2) = PopmeanCO(:,2);
% MSE(:,2)     = MSECO(:,2);

end