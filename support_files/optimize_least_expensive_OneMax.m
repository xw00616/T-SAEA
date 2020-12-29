function [Pop,Fitness] = optimize_least_expensive_OneMax(Population,Bounds,latency,Problem,corre,id_nex)
    
    FE_Max = size(Population,1)*(latency-1);
    
    Archive = call_GA(Population,FE_Max,Bounds,Problem,corre,id_nex);
    
    Pop = Archive(:,1:end-1);
    Fitness = Archive(:,end);
end

function Archive = call_GA(Population,FE_Max,Boundary,Problem,corre,id_nex)
    no_var = size(Boundary,2);
    F = P_objective1('value', Problem, corre, Population);
    FunctionValue = F(:,id_nex);
    FE = size(Population,1);
    Archive = [Population,FunctionValue];
    N = 100;
    while FE <FE_Max
        MatingPool = F_mating1(Population,N);
        Coding = 'Real';
        Offspring = P_generator(MatingPool,Boundary,Coding,N);         
        Offspring = unique(Offspring,'rows');
        Lia = ismember(Offspring,Archive(:,1:no_var),'rows');
        r_unique = find(Lia(:,1)==0);
        Offspring = Offspring(r_unique,:);        
        F_Values = P_objective1('value', Problem, corre, Offspring);
      
        FE = FE + size(Offspring,1);
        Fitness = F_Values(:,id_nex);
        Archive = [Archive;[Offspring,Fitness]];
        Population = [Population;Offspring];
        FunctionValue = [FunctionValue;Fitness];
        if(mod(size(FunctionValue,1),2) == 1)
            FunctionValue = [FunctionValue; FunctionValue(1,:)];
            Population = [Population;Population(1,:)];
        end
        selection = tournamentSelection(2,FunctionValue);
        Population = Population(selection(1,:),:);
        FunctionValue = FunctionValue(selection(1,:),:);     
   end
end