function P_output_OneMax (Population,time,Algorithm,Problem,corre,Run)

M=2;
FunctionValue = P_objective1('value',Problem,corre,Population);

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

eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue time',])
end


