function P_output(Population,Algorithm,Problem,M,Run)


FunctionValue = P_objective('value',Problem,M,Population);
if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end;
TrueValue = P_objective('true',Problem,M,1000);

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);
GDvalue = GD(FunctionValue,TrueValue);
Spreadvalue = Spread(FunctionValue,TrueValue);
[HVvalue,PopObj] = HV(FunctionValue,TrueValue);
IGDvalue = IGD(FunctionValue,TrueValue);
PF=size(FunctionValue,1);

if(M == 2)
    Plot2D(TrueValue, FunctionValue, 'ro');
end;



eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue HVvalue IGDvalue GDvalue Spreadvalue',])
end


