function [D,varargout] = P_settings(Algorithm,Problem,M)
    D = set_problem(Problem,M);
    Parameter = set_algorithm(Algorithm,Problem,M);
    varargout = num2cell(Parameter);
end

function D = set_problem(Problem,M)
    k = find(~isstrprop(Problem,'digit'),1,'last');%isstrprop??????????find?????????
    D1 = str2double(Problem(k+1:end));

    Problem = Problem(1:k);
    switch D1
        case {1,8}
            D=M+4;
        case {2,3,4,5,6,9}
            D=M+9;
        case 7
            D=M+19;
%         case {8,9}
%             D=10*M;
    end
                
    switch Problem
        case {'DTLZ', 'SDTLZ'}
            if M < 2 || M > 10
                error('Objective Number Not Supported !');
            end
            if D1 < 1 || D1 > 10
                error([Problem,'Not Exist']);
            end
         case {'UF'}
             D=30;
        case  {'ZDT'}
            D=8;
        case  {'WFG'}
            D=11;
            case  {'OneMax'}
            D=10;
        otherwise
            error([Problem,'Not Exist']);
    end
%                 Generations = [1000 500 1000 500 300 100 100 100 100 ];
%             Generations = Generations(D1);
%     N = [100 105 120 126 132 112 156 90 275];
%     N = N(M-1);
end

function Parameter = set_algorithm(Algorithm,Problem,M)
    Parameter = NaN;
    k = find(~isstrprop(Problem,'digit'),1,'last');
    D = str2double(Problem(k+1:end));
    Problem = Problem(1:k);
    switch Algorithm
        case {'RVEA'}
            p1 = [99 13  7  5  4  3  3  2  3];
            p2 = [ 0  0  0  0  1  2  2  2  2];
            p1 = p1(M-1);
            p2 = p2(M-1);
            Parameter(1) = p1; Parameter(2) = p2;
        case 'null'
        otherwise
            error([Algorithm,'Not Exist']);
    end
end

