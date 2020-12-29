function f  = evaluate_least_expensive_obj_OneMax(Population,Problem,id_nex)
    
    F = P_objective('value', Problem,2,Population);
    f = F(:,id_nex);

end