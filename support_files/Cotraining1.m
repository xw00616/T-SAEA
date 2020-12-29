function [Semi1Pop,Semi1Obj,Semi2Pop,Semi2Obj]=Cotraining1(PopNew1,PopNew1_TL,Model,FE)
k          = min(15,size(PopNew1,1));
PopSem     = [PopNew1;PopNew1_TL];
NumN       = size(PopSem,1);
Popmean    = zeros(NumN,2);
MSE        = zeros(NumN,2);
for j = 1 : 2
    [Popmean(:,j),MSE(:,j)] = predict(Model{j+1},PopSem);%   use mean
end

% a           = -0.5*cos(FE*pi/200)+0.5;
% b           = 1-a;
% [MMSE,~]    = max(MSE,[],1);
% [MPopObj,~] = max(Popmean,[],1);
% Popmeand    = [abs(Popmean(:,1)-Popmean(:,2)),abs(Popmean(:,2)-Popmean(:,1))];
% fit=Popmeand./repmat(MPopObj,size(Popmean,1),1)+MSE./repmat(MMSE,size(Popmean,1),1);
% fit         = Popmeand./repmat(MPopObj,size(Popmean,1),1)*b+MSE./repmat(MMSE,size(Popmean,1),1)*a;
% [fit1,I]    = sort(fit);
% Semi1Pop    = PopSem(I(1:k,2),:);
% Semi1Obj    = Popmean(I(1:k,2),2);
% Semi2Pop    = PopSem(I(1:k,1),:);
% Semi2Obj    = Popmean(I(1:k,1),1);

[SortUnc,I] = sort(MSE);
Semi1Pop    = PopSem(I(1:k,2),:);
Semi1Obj    = Popmean(I(1:k,2),2);
Semi2Pop    = PopSem(I(1:k,1),:);
Semi2Obj    = Popmean(I(1:k,1),1);
end