% The selection function in RVEA
function [Selection] = F_select(FunctionValue, V, theta0)

cosineVV = V*V';
[scosineVV, neighbor] = sort(cosineVV, 2, 'descend');%sort(A,2) sorts the elements of each row.
acosVV = acos(scosineVV(:,2));%Inverse cosine in radians
refV = (acosVV);

[N M] = size(FunctionValue);
VN = size(V, 1);

Zmin = min(FunctionValue,[],1);
Zmax = max(FunctionValue,[],1);

%Translation
FunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));

%Solutions associattion to reference vectors
clear class;
uFunctionValue = FunctionValue./repmat(sqrt(sum(FunctionValue.^2,2)), [1 M]);
cosine = uFunctionValue*V'; %calculate the cosine values between each solution and each vector
acosine = acos(cosine);
[maxc maxcidx] = max(cosine, [], 2);
class = struct('c', cell(1,VN)); %classification
for i = 1:N
    class(maxcidx(i)).c = [class(maxcidx(i)).c, i];
end;

Selection = [];
for k = 1:VN
    if(~isempty(class(k).c))
        sub = class(k).c;
        subFunctionValue = FunctionValue(sub,:);
        
        %APD calculation
        subacosine = acosine(sub, k);
        D1 = sqrt(sum(subFunctionValue.^2,2));
        subacosine = subacosine/refV(k);
        D = D1.*(1 + (theta0)*(subacosine));
        
        [mind mindidx] = min(D);
        Selection = [Selection; sub(mindidx)];
    end;
end;

end

