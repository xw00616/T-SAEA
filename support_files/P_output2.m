function P_output2 (Population,time,Algorithm,Problem,M,Run)

% M=2;
Population1=Population;
FunctionValue = P_objective1('value',Problem,M,Population);

% TrueValue = P_objective1('true',Problem,M,1000);

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);
% GDvalue = GD(FunctionValue,TrueValue);
% Spreadvalue = Spread(FunctionValue,TrueValue);
% HVvalue = HV(FunctionValue,TrueValue);
% IGDvalue = IGD(FunctionValue,TrueValue);
% PF=size(FunctionValue,1);

% if(M == 2)
%     Plot2D(TrueValue, FunctionValue, 'ro');
% end;

% if(M == 3)
%     Plot3D_PFboundry(TrueValue, FunctionValue, 'ro');
% end;

eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(Run),' Population1 Population FunctionValue time',])
end


