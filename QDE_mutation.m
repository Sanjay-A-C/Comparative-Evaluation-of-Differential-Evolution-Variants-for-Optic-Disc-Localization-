function dv = QDE_mutation(pop, F, r1, r2, r3)
    
    

    % Quantum Perturbation using rotation-based update
    theta = rand() * pi; % Quantum rotation angle
     %theta = rand(1, 4) * 2 * pi;
    Q_theta = sin(theta).^2; % Quantum perturbation factor

    % Differential Evolution Mutation with Quantum Influence
    %dv = round(abs(pop(r1,:)+F*(-(pop(r2,:)-pop(r3,:)))+ Q_theta * randn(size(pop(r1,:))))) ;

    dv = round(abs(pop(r1,:)+F*(-(pop(r2,:)-pop(r3,:)))+ Q_theta)) ;
   
    
    % Check and correct out-of-bound values:
    if dv(1) < 150 || dv(1) > 1350
        dv(1) = pop(r1, 1); 
    end
    if dv(2) < 150 || dv(2) > 1002
        dv(2) = pop(r1, 2);
    end
    if dv(3) < -10 || dv(3) > 10
        dv(3) = pop(r1, 3); 
    end
    if dv(4) < -10 || dv(4) > 10
        dv(4) = pop(r1, 4); 
    end
end



