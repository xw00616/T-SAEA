clear;clc;
itr=20;Ke=5;
Problems={'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'};
strnum=2;V=20;iteration=24;M=3;

for Prob = 1:length(Problems)
    clear PF;
    Problem=Problems{Prob};
    array=[];
    k = find(~isstrprop(Problem,'digit'),1,'last');%k是Problem中字母位数
    str1=[Problem(1) Problem(k+1) 'KRVEA'];
    str2=[Problem(1) Problem(k+1) 'KRVEA1'];
    %=[Problem(1) Problem(k+1) 'KRVEA2'];
    %str3=[Problem(1) Problem(k+1) 'KRVEA42'];
    %str3=[Problem(1) Problem(k+1) 'KRVEA41'];
    %str4=[Problem(1) Problem(k+1) 'KRVEA3'];
    load(str1);load(str2);%load(str3);%load(str4);%load(str5);load(str6);load(str7);load(str8);load(str9);
    %% 求参考点
    num=1000;
    PF=WFG_PF(num,Problem,M);
    %% 画HV随代数变化图
    N=zeros(strnum,iteration+1);
    for j=1:strnum
        NN=zeros(itr,1);
        a=eval(eval(['str',num2str(j)]));
        for k=1:itr
            aa=a(k).ch;aa=aa(:,V+1:V+M);
            for i=1:iteration+1
                b=aa(1:(11*V-1+Ke*(i-1)),:);
                b=non_sort(b,M);bb=b(find(b(:,M+1)==1),1:M);
                igd=P_evaluate('IGD',bb,PF);NN(k,i)=igd;
            end
        end
        r(j).NN=NN;
        N(j,:)=mean(NN,1);
        BZC(j,:)=std(NN(:,end));
        disp(sprintf('Prob %u : %u / %u method finished', Prob, j, strnum));
    end
    a=['Prob',num2str(Prob),'N','.mat'];save(a,'N')    
    a=['Prob',num2str(Prob),'NN','.mat'];save(a,'r')   
    a=['Prob',num2str(Prob),'std','.mat'];save(a,'BZC')
end
