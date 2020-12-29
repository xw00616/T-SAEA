function [PopDec,PopObj]=TrainingSelect(PopDec,PopObj,L,NumberNew)
                PopDec1  =  PopDec(1:NumberNew,:);
                PopObj1  =  PopObj(1:NumberNew,:);
                [~,distinct1]  = unique(PopDec,'rows');
                PopDec    = PopDec(distinct1,:);
                PopObj    = PopObj(distinct1,1);
                datalong=size(PopDec,1);
            if datalong <=L
                disp(sprintf('No training data decrease'));
                PopDec=PopDec;PopObj=PopObj;
            else
                
                [val,paixu] = sort(PopObj);
                data=paixu(1:floor(L/2), :);
                paixu=paixu(floor(L/2)+1:end,:);
                index=randperm(size(paixu,1));
                data=[data;paixu(index(1:L-floor(L/2)-NumberNew), :)];
                PopDec=[PopDec(data,:);PopDec1];PopObj=[PopObj(data,:);PopObj1];
            end
    end
