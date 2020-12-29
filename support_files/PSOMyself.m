function index=PSOMyself(x, y)

global D
m=20;%种群数量
iter=30;%迭代次数
threshold=0.5;%CSO参数
c1=2;c2=2;w=1.05;%CSO参数
%position上下限
lu = [zeros(1, D); ones(1, D)];
XRRmin = repmat(lu(1, :), m, 1);
XRRmax = repmat(lu(2, :), m, 1);
VelMax=(XRRmax-XRRmin)/10;

%position
p = rand(m, D);
p=XRRmin + (XRRmax - XRRmin) .*p;
p(find(p<threshold))=0;p(find(p>=threshold))=1;
v = VelMax .* rand(m,D);

fitness = PSO_CostFunction(p, x, y);

[bestever, index] = min(fitness);
zbest=p(index,:);fitnesszbest=bestever;
gbest=p;fitnessgbest=fitness;

for i = 1 : iter
    
    v = w*v + c1*rand(m,D).*(gbest-p) + c2*rand(m,D).*(repmat(zbest, m, 1)-p);
    v = max(min(v,VelMax), -VelMax);
    p = p + v;
    p(find(p<threshold))=0;p(find(p>=threshold))=1;
   
    fitness = PSO_CostFunction(p, x, y); 
    
    index = find(fitness < fitnessgbest);
    fitnessgbest(index) = fitness(index);
    gbest(index,:) = p(index,:);
    [bestever, index] = min(fitness);
    if  bestever < fitnesszbest
        fitnesszbest = bestever;
        zbest = p(index, :);
    end
end;
index=find(zbest==1);





    

