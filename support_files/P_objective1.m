%This file includes the original and scaled DTLZ1 to DTLZ4

function [Output,Boundary,Coding] = P_objective1(Operation,Problem,M,Input)
k = find(~isstrprop(Problem,'digit'),1,'last');
switch Problem(1:k)
    case 'DTLZ'
        [Output,Boundary,Coding] = P_DTLZ(Operation,Problem,M,Input);
    case 'UF'
        [Output,Boundary,Coding] = P_UF(Operation,Problem,M,Input);
    case 'ZDT'
        [Output,Boundary,Coding] = P_ZDT(Operation,Problem,M,Input);
    case 'WFG'
        [Output,Boundary,Coding] = P_WFG(Operation,Problem,M,Input);
    case 'OneMax'
        [Output,Boundary,Coding] =P_OneMax(Operation,Problem,M,Input);
    otherwise
        error([Problem,'Not Exist']);
end
end

function [Output,Boundary,Coding] = P_DTLZ(Operation,Problem,M,Input)
persistent K;
Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    case 'init'
        k = find(~isstrprop(Problem,'digit'),1,'last');
        K = [5 10 10 10 10 10 20 5 10];
        K = K(str2double(Problem(k+1:end)));
        
        D = M+K-1;
        
        MaxValue   = ones(1,D);
        MinValue   = zeros(1,D);
        %         Population = rand(Input,D);
        %         Population = lhsamp(Input,D);
        %             Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
        Population= repmat(MinValue,Input,1) + repmat(MaxValue-MinValue,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
        
        Output   = Population;
        Boundary = [MaxValue;MinValue];
        Coding   = 'Real';
        %Objective Function Evaluation
    case 'value'
        k = find(~isstrprop(Problem,'digit'),1,'last');
        K = [5 10 10 10 10 10 20 5 10];
        K = K(str2double(Problem(k+1:end)));
        Population    = Input;
        FunctionValue = zeros(size(Population,1),M);
        switch Problem
            case 'DTLZ8'      %DTLZ1a
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(2.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                    end
                end
            case 'DTLZ1'
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                    end
                end
            case 'DTLZ2'
                g = sum((Population(:,M:end)-0.5).^2,2);
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'DTLZ3'
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'DTLZ9'%DTLZ3a
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(2.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'DTLZ4'
                Population(:,1:M-1) = Population(:,1:M-1).^100;
                g = sum((Population(:,M:end)-0.5).^2,2);
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'DTLZ5'
                g      = sum((Population(:,M:end)-0.5).^2,2);
                Temp   = repmat(g,1,M-2);
                Population(:,2:M-1) = (1+2*Temp.*Population(:,2:M-1))./(2+2*Temp);
                FunctionValue = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(Population(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(Population(:,M-1:-1:1)*pi/2)];
            case 'DTLZ6'
                
                g      = sum(Population(:,M:end).^0.1,2);
                Temp   = repmat(g,1,M-2);
                Population(:,2:M-1) = (1+2*Temp.*Population(:,2:M-1))./(2+2*Temp);
                FunctionValue = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(Population(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(Population(:,M-1:-1:1)*pi/2)];
            case 'DTLZ7'
                FunctionValue         = zeros(size(Population,1),M);
                g               = 1+9*mean(Population(:,M:end),2);
                FunctionValue(:,1:M-1) = Population(:,1:M-1);
                FunctionValue(:,M)     = (1+g).*(M-sum(Population(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*FunctionValue(:,1:M-1))),2));
                %             case 'DTLZ10'
                %                 g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                %                 FunctionValue(:,1) = 0.5.*prod(Population(:,1:M-1),2).*(1+g);
                %                 FunctionValue(:,2) = 0.5.*prod(Population(:,1:M-1),2).*(1+g);
                %                     FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                %
                %                 g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(2.*pi.*(Population(:,M:end)-0.5)),2));
                %                 for i = 1 : M
                %                     FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                %                     if i > 1
                %                         FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                %                     end
                %                 end
                %                 %             case 'DTLZ9'
                %                 [N,D]  = size(Population);
                %                 Population = Population.^0.1;
                %
                %                 FunctionValue = zeros(N,M);
                %                 for m = 1 : M
                %                     FunctionValue(:,m) = sum(Population(:,(m-1)*D/M+1:m*D/M),2);
                %                 end
        end
        Output = FunctionValue;
        %Sample True PFs
    case 'true'
        switch Problem
            case {'DTLZ1','DTLZ8'}
                Population = T_uniform(Input,M);
                Population = Population/2;
            case {'DTLZ2','DTLZ3','DTLZ4','DTLZ9'}
                Population = T_uniform(Input,M);
                %                 for i = 1 : size(Population,1)
                %                     Population(i,:) = Population(i,:)./norm(Population(i,:));
                %                 end
                Population =Population ./repmat(sqrt(sum(Population .^2,2)),1,M);
            case {'DTLZ5','DTLZ6'}
                N=Input;
                Population = [0:1/(N-1):1;1:-1/(N-1):0]';
                Population = Population./repmat(sqrt(sum(Population.^2,2)),1,size(Population,2));
                Population= [Population(:,ones(1,M-2)),Population];
                Population= Population./sqrt(2).^repmat([M-2,M-2:-1:0],size(Population,1),1);
            case 'DTLZ7'
                N=Input;
                interval     = [0,0.251412,0.631627,0.859401];
                median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
                X            = ReplicatePoint(N,M-1);
                X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
                X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
                Population        = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
                
        end
        Output = Population;
end
end
function W = ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end


function W = T_uniform(k,M)
H = floor((k*prod(1:M-1))^(1/(M-1)));
while nchoosek(H+M-1,M-1) >= k && H > 0
    H = H-1;
end
if nchoosek(H+M,M-1) <= 2*k || H == 0
    H = H+1;
end
k = nchoosek(H+M-1,M-1);
Temp = nchoosek(1:H+M-1,M-1)-repmat(0:M-2,nchoosek(H+M-1,M-1),1)-1;
W = zeros(k,M);
W(:,1) = Temp(:,1)-0;
for i = 2 : M-1
    W(:,i) = Temp(:,i)-Temp(:,i-1);
end
W(:,end) = H-Temp(:,end);
W = W/H;
end


function [Output,Boundary,Coding] = P_UF(Operation,Problem,M,Input)

Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    case 'init'
        D = 30;
        N=Input;
        switch Problem
            case  {'UF1','UF2','UF5','UF6','UF7'}
                lower    = [0,zeros(1,D-1)-1];
                upper    = ones(1,D);
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
            case 'UF3'
                lower    = zeros(1,D);
                upper    = ones(1,D);
                %                  Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
            case 'UF4'
                lower    = [0,zeros(1,D-1)-2];
                upper    = [1,zeros(1,D-1)+2];
                %                  Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
                %Objective Function Evaluation
        end
    case 'value'
        X   = Input;
        PopObj = zeros(size(X,1),M);
        switch Problem
            case 'UF1'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
                PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
                
            case 'UF2'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y       = zeros(size(X));
                X1      = repmat(X(:,1),1,length(J1));
                Y(:,J1) = X(:,J1)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
                X1      = repmat(X(:,1),1,length(J2));
                Y(:,J2) = X(:,J2)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D);
                PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
                PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
            case 'UF3'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - repmat(X(:,1),1,D).^(0.5*(1+3*(repmat(1:D,size(X,1),1)-2)/(D-2)));
                PopObj(:,1) = X(:,1)         + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
                PopObj(:,2) = 1-sqrt(X(:,1)) + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
            case 'UF4'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                hY = abs(Y)./(1+exp(2*abs(Y)));
                PopObj(:,1) = X(:,1)      + 2*mean(hY(:,J1),2);
                PopObj(:,2) = 1-X(:,1).^2 + 2*mean(hY(:,J2),2);
            case 'UF5'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                hY = 2*Y.^2 - cos(4*pi*Y) + 1;
                PopObj(:,1) = X(:,1)   + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J1),2);
                PopObj(:,2) = 1-X(:,1) + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J2),2);
            case 'UF6'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                PopObj(:,1) = X(:,1)   + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
                PopObj(:,2) = 1-X(:,1) + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
                
            case 'UF7'
                D  = size(X,2);
                J1 = 3 : 2 : D;
                J2 = 2 : 2 : D;
                Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                PopObj(:,1) = X(:,1).^0.2   + 2*mean(Y(:,J1).^2,2);
                PopObj(:,2) = 1-X(:,1).^0.2 + 2*mean(Y(:,J2).^2,2);
                
        end
        Output =  PopObj;
        %Sample True PFs
    case 'true'
        N=Input;
        switch Problem
            case {'UF1','UF2','UF3'}
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1).^0.5;
            case 'UF4'
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1).^2;
            case 'UF5'
                P(:,1) = (0:1:20)'/20;
                P(:,2) = 1 - P(:,1);
            case 'UF6'
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1);
                P(P(:,1)>0 & P(:,1)<1/4 | P(:,1)>1/2 & P(:,1)<3/4,:) = [];
            case 'UF7'
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1);
        end
        Output = P;
        
end
end

function [Output,Boundary,Coding] = P_ZDT(Operation,Problem,M,Input)

Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    case 'init'
        D = 8;
        N=Input;
        switch Problem
            case  {'ZDT1','ZDT2','ZDT3','ZDT6'}
                %                 obj.Global.M = 2;
                lower    = zeros(1,D);
                upper    = ones(1,D);
                %                 Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
            case 'ZDT4'
                lower    = [0,zeros(1,D-1)-5];
                upper    = [1,zeros(1,D-1)+5];
                %                 Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
        end
        
        %Objective Function Evaluation
    case 'value'
        PopDec   = Input;
        PopObj = zeros(size(PopDec,1),M);
        
        switch Problem
            case 'ZDT1'
                PopObj(:,1) = PopDec(:,1);
                g = 1 + 9*sum(PopDec(:,2:end),2);
                h = 1 - (PopObj(:,1)./g).^0.5;
                PopObj(:,2) = g.*h;
            case 'ZDT2'
                PopObj(:,1) = PopDec(:,1);
                g = 1 + 9*sum(PopDec(:,2:end),2);
                h = 1 - (PopObj(:,1)./g).^2;
                PopObj(:,2) = g.*h;
            case 'ZDT3'
                PopObj(:,1) = PopDec(:,1);
                g = 1 + 9*sum(PopDec(:,2:end),2);
                h = 1 - (PopObj(:,1)./g).^0.5 - PopObj(:,1)./g.*sin(10*pi*PopObj(:,1));
                PopObj(:,2) = g.*h;
            case 'ZDT4'
                PopObj(:,1) = PopDec(:,1);
                g = 1 + 10*(size(PopDec,2)-1) + sum(PopDec(:,2:end).^2-10*cos(4*pi*PopDec(:,2:end)),2);
                h = 1 - (PopObj(:,1)./g).^0.5;
                PopObj(:,2) = g.*h;
            case 'ZDT6'
                PopObj(:,1) = 1 - exp(-4*PopDec(:,1)).*sin(6*pi*PopDec(:,1)).^6;
                g = 1 + 9*sum(PopDec(:,2:end),2).^0.25;
                h = 1 - (PopObj(:,1)./g).^2;
                PopObj(:,2) = g.*h;
        end
        Output = PopObj;
        %Sample True PFs
    case 'true'
        N=Input;
        switch Problem
            case {'ZDT1','ZDT4'}
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1).^0.5;
            case 'ZDT2'
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1).^2;
            case 'ZDT3'
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = 1 - P(:,1).^0.5 - P(:,1).*sin(10*pi*P(:,1));
                P      = P(NDSort(P,1)==1,:);
            case 'ZDT6'
                minf1  = 0.280775;
                P(:,1) = (minf1:(1-minf1)/(N-1):1)';
                P(:,2) = 1 - P(:,1).^2;
        end
        Output = P;
end
end
function [Output,Boundary,Coding] = P_WFG(Operation,Problem,M,Input)
persistent K;
K=1;
D = M+9;
Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    
    case 'init'
        N=Input;
        switch Problem
            case  {'WFG1'}
                lower    = zeros(1,D);
                upper    = 2 : 2 : 2*D;
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
            case {'WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'}
                %                 D        = ceil((D-K)/2)*2 + K;
                lower    = zeros(1,D);
                upper    = 2 : 2 : 2*D;
                Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [upper;lower];
                Coding   = 'Real';
        end
        
        %Objective Function Evaluation
    case 'value'
        PopDec   = Input;
        PopObj = zeros(size(PopDec,1),M);
        
        switch Problem
            case 'WFG1'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1(:,1:K)     = z01(:,1:K);
                t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
                
                t2 = zeros(N,K+L);
                t2(:,1:K)     = t1(:,1:K);
                t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);
                
                t3 = zeros(N,K+L);
                t3 = b_poly(t2,0.02);
                
                t4 = zeros(N,M);
                for i = 1 : M-1
                    t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
                end
                t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
                end
                x(:,M) = t4(:,M);
                
                h      = convex(x);
                h(:,M) = mixed(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WFG2'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1(:,1:K)     = z01(:,1:K);
                t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
                
                t2 = zeros(N,K+L/2);
                t2(:,1:K) = t1(:,1:K);
                % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
                t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
                % ---------------------------------------------------------
                
                t3 = zeros(N,M);
                for i = 1 : M-1
                    t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
                end
                x(:,M) = t3(:,M);
                
                h      = convex(x);
                h(:,M) = disc(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WFG3'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = [1,zeros(1,M-2)];
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1(:,1:K)     = z01(:,1:K);
                t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
                
                t2 = zeros(N,K+L/2);
                t2(:,1:K) = t1(:,1:K);
                % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
                t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
                % ---------------------------------------------------------
                
                t3 = zeros(N,M);
                for i = 1 : M-1
                    t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
                end
                x(:,M) = t3(:,M);
                
                h      = linear(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WFG4'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1 = s_multi(z01,30,10,0.35);
                
                t2 = zeros(N,M);
                for i = 1 : M-1
                    t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
                end
                x(:,M) = t2(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WGF5'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1 = s_decept(z01,0.35,0.001,0.05);
                
                t2 = zeros(N,M);
                for i = 1 : M-1
                    t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
                end
                x(:,M) = t2(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WGF6'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1(:,1:K)     = z01(:,1:K);
                t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
                
                t2 = zeros(N,M);
                for i = 1 : M-1
                    t2(:,i) = r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
                end
                % Same as <t2(:,M)=r_nonsep(t1(:,K+1:end),L)>
                SUM = zeros(N,1);
                for i = K+1 : K+L-1
                    for j = i+1 : K+L
                        SUM = SUM + abs(t1(:,i)-t1(:,j));
                    end
                end
                t2(:,M) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
                % -------------------------------------------
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
                end
                x(:,M) = t2(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WGF7'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
                Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
                t1(:,1:K) = z01(:,1:K).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K)).*abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
                % ------------------------------------------------------------------------------------------
                t1(:,K+1:end) = z01(:,K+1:end);
                
                t2 = zeros(N,K+L);
                t2(:,1:K)     = t1(:,1:K);
                t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);
                
                t3 = zeros(N,M);
                for i = 1 : M-1
                    t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
                end
                x(:,M) = t3(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WGF8'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                t1(:,1:K) = z01(:,1:K);
                % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,1:i-1),ones(1,i-1)),0.98/49.98,0.02,50)>
                Y = (cumsum(z01,2)-z01)./repmat(0:K+L-1,N,1);
                t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,K+1:K+L)).*abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98)));
                % --------------------------------------------------------------------------------------
                
                t2 = zeros(N,K+L);
                t2(:,1:K)     = t1(:,1:K);
                t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);
                
                t3 = zeros(N,M);
                for i = 1 : M-1
                    t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
                end
                t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
                end
                x(:,M) = t3(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            case 'WGF9'
                [N,D] = size(PopDec);
                L = D - K;
                D = 1;
                S = 2 : 2 : 2*M;
                A = ones(1,M-1);
                
                z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
                
                t1 = zeros(N,K+L);
                % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
                Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
                t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
                % ------------------------------------------------------------------------------------------
                t1(:,end)     = z01(:,end);
                
                t2 = zeros(N,K+L);
                t2(:,1:K)     = s_decept(t1(:,1:K),0.35,0.001,0.05);
                t2(:,K+1:end) = s_multi(t1(:,K+1:end),30,95,0.35);
                
                t3 = zeros(N,M);
                for i = 1 : M-1
                    t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
                end
                % Same as <t3(:,M)=r_nonsep(t2(:,K+1:end),L)>
                SUM = zeros(N,1);
                for i = K+1 : K+L-1
                    for j = i+1 : K+L
                        SUM = SUM + abs(t2(:,i)-t2(:,j));
                    end
                end
                t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
                % -------------------------------------------
                
                x = zeros(N,M);
                for i = 1 : M-1
                    x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
                end
                x(:,M) = t3(:,M);
                
                h = concave(x);
                PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
        end
        Output = PopObj;
        %Sample True PFs
    case 'true'
        N=Input;
        switch Problem
            case {'WFG1'}
                P = UniformPoint(N,M);
                c = ones(size(P,1),M);
                for i = 1 : size(P,1)
                    for j = 2 : M
                        temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                        c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                    end
                end
                x = acos(c)*2/pi;
                temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
                a = 0 : 0.0001 : 1;
                E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
                [~,rank] = sort(E,2);
                for i = 1 : size(x,1)
                    x(i,1) = a(min(rank(i,1:10)));
                end
                P      = convex(x);
                P(:,M) = mixed(x);
                P      = repmat(2:2:2*M,size(P,1),1).*P;
            case 'WFG2'
                P = UniformPoint(N,M);
                c = ones(size(P,1),M);
                for i = 1 : size(P,1)
                    for j = 2 : M
                        temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                        c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                    end
                end
                x = acos(c)*2/pi;
                temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
                a = 0 : 0.0001 : 1;
                E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
                [~,rank] = sort(E,2);
                for i = 1 : size(x,1)
                    x(i,1) = a(min(rank(i,1:10)));
                end
                P      = convex(x);
                P(:,M) = disc(x);
                P      = P(NDSort(P,1)==1,:);
                P      = repmat(2:2:2*M,size(P,1),1).*P;
            case 'WFG3'
                X = (0:1/(N-1):1)';
                X = [X,zeros(N,M-2)+0.5,zeros(N,1)];
                P = linear(X);
                P = repmat(2:2:2*M,size(P,1),1).*P;
            case {'WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'}
                P = UniformPoint(N,M);
                P = P./repmat(sqrt(sum(P.^2,2)),1,M);
                P = repmat(2:2:2*M,size(P,1),1).*P;
        end
        Output = P;
end
end
%%%%1
function Output = s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
Output = roundn(Output,-6);
end

function Output = b_poly(y,a)
Output = y.^a;
end

function Output = r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end
function Output = mixed(x)
Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end
%%%%2
function Output = r_nonsep(y,A)
Output = zeros(size(y,1),1);
for j = 1 : size(y,2)
    Temp = zeros(size(y,1),1);
    for k = 0 : A-2
        Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
    end
    Output = Output+y(:,j)+Temp;
end
Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end
function Output = disc(x)
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end
%%3
function Output = linear(x)
Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end
%%4
function Output = s_multi(y,A,B,C)
Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
%%5
function Output = s_decept(y,A,B,C)
Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end
%%6

function [Output,Boundary,Coding] = P_OneMax(Operation,Problem,corre,Input)
persistent K;
K=1;
D =10;
M=2;
Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    
    case 'init'
        N=Input;
        lower    = zeros(1,D);
        upper    = ones(1,D);
        Population= repmat( lower,Input,1) + repmat(upper -lower,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
        Output   = Population;
        Boundary = [upper;lower];
        Coding   = 'Real';
        
        
        %Objective Function Evaluation
    case 'value'
        PopDec   = Input;
        mapx     = Input;
        PopObj = zeros(size(PopDec,1),M);
        % sum(x)
        PopObj(:,1)=sum(PopDec,2);
        % sum(y)
        for i=1:size(Input,1)
            for j=1:size(Input,2)
                pro=0.5*(1+corre);
                if rand<pro
                    mapx(i,j)=0;
                else
                    mapx(i,j)=1;
                end
            end
        end
        PopObj(:,2)=sum(abs(Input-mapx),2);
        Output = -PopObj;
        %Sample True PFs
    case 'true'
        Output=Input;
        
end
end




