function f  = evaluate_least_expensive_obj(Population,Problem,id_nex)
    
    F = P_objective('value', Problem,M,Population);
    f = F(:,id_nex);

end