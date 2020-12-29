function [PopNew,PopNewObj,PopNew1,PopNewObj1]=IndividualSelect(PopDec, PopObj, MSE,FE,V1,NewNum)

            a=-0.5*cos(FE*pi/200)+0.5;
            b=1-a;
            [MMSE,~]=max(MSE,[],1);
            [MPopObj,~]=max(PopObj,[],1);                 
            fit=PopObj./repmat(MPopObj,size(PopObj,1),1)*b+MSE./repmat(MMSE,size(PopObj,1),1)*a;
            %%%%%%%%%%select by the reference vectors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Selection = FSelection(fit,V1,(FE/200)^2);
            PopNewSelected = PopDec(Selection,:);
            PopObjSelected= PopObj(Selection,:);
            if size(PopNewSelected,1)<NewNum*5
                PopNew1=[PopNewSelected];
                PopNewObj1=[PopObjSelected];
            elseif size(PopNewSelected,1)>=NewNum*5
                aa=randperm(size(PopNewSelected,1),4);
                PopNew1=PopNewSelected(aa,:);
                PopNewObj1=PopObjSelected(aa,:);
            end
            Nondominate = P_sort(PopObjSelected,'first')==1;
            PopNewSelected=PopNewSelected(Nondominate,:);
            PopObjSelected=PopObjSelected(Nondominate,:);
            if size(PopNewSelected,1)<NewNum
                PopNew=PopNewSelected;
                PopNewObj=PopObjSelected;
            elseif size(PopNewSelected,1)>=NewNum
                aa=randperm(size(PopNewSelected,1),NewNum);
                PopNew=PopNewSelected(aa,:);
                PopNewObj=PopObjSelected(aa,:);
            end
end