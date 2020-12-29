%%GECCO paper%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%feature selection
%%%%%%%%%%%%%updating 3 per evaluation -------parameter(final version)
%%%%%%%%%%%%%model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all
addpath RVEA;
addpath Public;
global D
Algorithm = 'RVEA';
Pop_num=100;
objective=2;
Evaluations=200;
RunNum =20;
alpha=2;
Ratio=5;


%%%?????%%%%%%%%%%%%
Problems={'ZDT4','DTLZ8','DTLZ1','DTLZ4','DTLZ2','DTLZ3','DTLZ5','DTLZ6','DTLZ7','DTLZ9'...
    'UF1','UF2','UF3','UF4','UF5','UF6','UF7',...
    'ZDT1','ZDT2','ZDT3','ZDT4','ZDT6'};
% for Prob = 1:length(Problems)
for Prob =2
    Problem=Problems{Prob}
    tic
    for Run=5
        Run
        %% Generate random population
        M=objective;
        ArchiveAll=[];
        [D,p1,p2] = P_settings('RVEA',Problem,M);
%         upboundry=ones(1,D);
%         lpboundry=zeros(1,D);
%         Boundary=[upboundry;lpboundry];
        L    = 200;
        dd=ceil(0.5*D);
        N=100;
        THETA = 5.*ones(M,D);
        Model = cell(1,M);
        [Population,Boundary,Coding] = P_objective1('init',Problem,M,N);
         upboundry=Boundary(1,:);
        lpboundry=Boundary(2,:);
        FunctionValue = P_objective1('value',Problem,M,Population);
        selection=PSOMyself(Population, FunctionValue);
        FunctionValue2 = FunctionValue(:,2);
        Population2=Population;
        [V0,Pop_num] = UniformPoint(Pop_num,M);
        V     = V0;
        V10    = [V0];
        FE = size(Population,1);
        wmax=20;
        iter=1;
        Best_pos=0.5*ones(1,D);
        PopNew=[];
        Achx2=[];
        Achf2=[];
        MeanPredict=[];
        MeanTrue=[];
        MeanValue=[];
        MeanNewValue=[];
        %% Optimization
        while FE<=Evaluations
            %%train GP for expensive one
            %                 Delete duplicated solutions
             a=-0.5*cos(FE*pi/200)+0.5;
             b=1-a;
            if (mod(iter, Ratio) == 0||iter==1)
                [~,distinct1]  = unique(Population,'rows');
                PopDec2 = Population(distinct1,:);
                PopObj2 = FunctionValue(distinct1,2);
                datalong=size(PopDec2,1);
                if datalong <=L
                    disp(sprintf('No training data decrease'));
                    PopDec2=PopDec2;PopObj2=PopObj2;
                else
                    [val2,paixu2] = sort(PopObj2);
                    data2=paixu2(1:floor(L/2), :);
                    paixu2=paixu2(floor(L/2)+1:end,:);
                    index=randperm(size(paixu2,1));
                    data2=[data2;paixu2(index(1:L-floor(L/2)), :)];
                    PopDec2=PopDec2(data2,:);PopObj2=PopObj2(data2,:);
                end
                
                dmodel     = dacefit(PopDec2,PopObj2,'regpoly0','corrgauss',THETA(2,:),1e-5.*ones(1,D),100.*ones(1,D));
                Model{2}   = dmodel;
                THETA(2,:) = dmodel.theta;
            else
                selection=PSOMyself(Population, FunctionValue);%%%%2 0.04 1 23
                Model{2}.theta(selection)=a.*THETA(2,selection)+b.*THETA(1,selection);%%%2 0.04 1 46
                THETA(2,selection)=a.*THETA(2,selection)+b.*THETA(1,selection);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%cheap one evaluate more solutions
            
            if iter==1
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [X_nex,F_nex] = optimize_least_expensive(Population,Boundary,Ratio,Problem,M,1);
                 Population1=X_nex;
                ArchiveAll=X_nex;
                FunctionValue1=F_nex;

                %%%%%train GP for the first time
                [~,index]  = unique(Population1,'rows');
                PopDec1 = Population1(index,:);
                PopObj1 = FunctionValue1(index,1);
                dmodel     = dacefit(PopDec1,PopObj1,'regpoly1','corrgauss',THETA(1,:),1e-5.*ones(1,D),100.*ones(1,D));
                Model{1}   = dmodel;
                THETA(1,:) = dmodel.theta;
            else
                [~,index]  = unique(Population1,'rows');
                PopDec1 = Population1(index,:);
                PopObj1 = FunctionValue1(index,1);
                datalong=size(PopDec1,1);
                if datalong <=L
                    disp(sprintf('No training data decrease'));
                    PopDec1=PopDec1;PopObj1=PopObj1 ;
                else
                    [val1,paixu1] = sort(PopObj1);
                    data1=paixu1(1:floor(L/2), :);
                    paixu1=paixu1(floor(L/2)+1:end,:);
                    index=randperm(size(paixu1,1));
                    data1=[data1;paixu1(index(1:L-floor(L/2)), :)];
                    PopDec1=PopDec1(data1,:);PopObj1=PopObj1(data1,:);
                end
                dmodel     = dacefit(PopDec1,PopObj1,'regpoly1','corrgauss',THETA(1,:),1e-5.*ones(1,D),100.*ones(1,D));
                Model{1}   = dmodel;
                THETA(1,:) = dmodel.theta;
            end
        
            %%the predict value
            %%RVEA optimization
            [V0,Pop_num] = UniformPoint(Pop_num,M);
            V10    = [V0];
            w=1;
            V1    = V10;
            cal=0;
            Pop=Population;
            while w < 20
                if w==1
                    Popmean = zeros(size(ArchiveAll,1),M);
                    MSE    = zeros(size(ArchiveAll,1),M);
                    for i = 1: size(ArchiveAll,1)
                        for j = 1 : M
                            [Popmean(i,j),~,MSE(i,j)] = predictor(ArchiveAll(i,:),Model{j});%   use mean
                        end
                    end
                    PopObj=Popmean;
                    Selection = FSelection(PopObj,V,(w/wmax)^alpha);
                    Pop = ArchiveAll(Selection,:);
                    Popmean=Popmean(Selection,:);
                    PopObj = PopObj(Selection,:);
                    MSE=MSE(Selection,:);
                else
                      [MatingPool] = F_mating(Pop,N);
                      OffDec = P_generator(MatingPool,Boundary,'Real',N);

                    Pop = [Pop;OffDec];
                    [N1,~]  = size(Pop);
                    Popmean = zeros(N1,M);
                    MSE    = zeros(N1,M);
                    for i = 1: N1
                        for j = 1 : M
                            [Popmean(i,j),~,MSE(i,j)] = predictor(Pop(i,:),Model{j});%   use mean
                        end
                    end
                    PopObj=Popmean;
                    Selection = FSelection(PopObj,V,(w/wmax)^alpha);
                    Pop = Pop(Selection,:);
                    Popmean=Popmean(Selection,:);
                    PopObj = PopObj(Selection,:);
                    MSE=MSE(Selection,:);
                end
                if(mod(w, ceil(wmax*0.1)) == 0)
                    V = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
                end;
                cal=cal+size(Pop,1);
                w=w+1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            a=-0.5*cos(FE*pi/200)+0.5;
            b=1-a;
            [MMSE,~]=max(MSE,[],1);
            [MPopObj,~]=max(PopObj,[],1);
            Nondominate = P_sort(FunctionValue,'first')==1;
            Best_pos=Population(Nondominate,:);
            fit=PopObj./repmat(MPopObj,size(PopObj,1),1)*b+MSE./repmat(MMSE,size(PopObj,1),1)*a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%original APD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Selection = FSelection(fit,V1,(FE/200)^alpha);
            PopNewSelected = Pop(Selection,:);
            if size(PopNewSelected,1)<3
                PopNew0=PopNewSelected(:,1:D);
            elseif size(PopNewSelected,1)>=3
                PopNew0=PopNewSelected(randperm(size(PopNewSelected,1),3),1:D);
            end
            if(mod(FE, ceil(200*0.1)) == 0)
                V1 = V10.*repmat(max(FunctionValue,[],1)-min(FunctionValue,[],1),size(V10,1),1);
            end;       
            %%%%%%%%%%%%%%%%sample by LHS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N1=3*Ratio;
            MinValue=min(PopNew0);
            MaxValue=max(PopNew0);
%             PopNew1= repmat(MinValue,N1,1) + repmat(MaxValue-MinValue,N1,1).*lhsdesign(N1,D,'criterion','maximin','iterations',1000);
           
            PopNew1   = (repmat(MaxValue-MinValue,N1,1).*lhsamp(N1,D)+repmat(MinValue,N1,1));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PopNew=[PopNew;PopNew0];
            PopNew2=[PopNew0;PopNew1];
            AddFunctionValue0=P_objective1('value',Problem,M,PopNew0);
            AddFunctionValue1=P_objective1('value',Problem,M,PopNew1);
            AddFunctionValue2=[AddFunctionValue0;AddFunctionValue1];

            
            Population1=[Population1;PopNew2];
            FunctionValue1=[FunctionValue1;AddFunctionValue2(:,1)];
            ArchiveAll=[ArchiveAll;PopNew2];
            
            
            FunctionValue=[FunctionValue;AddFunctionValue0];
            Population=[Population;PopNew0];
            FE = FE+size(PopNew0,1)
            P_output2(Population,toc,'TSAEA1000',Problem,M,Run);
            iter=iter+1;
        end
        Run=Run+1;
    end
end

