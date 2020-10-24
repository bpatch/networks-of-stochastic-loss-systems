function value = ExactEvaluationModelOne(C1, C2, Weights, Mu, lambda, CapacityCost)
%     C1 = 20;
%     C2 = 24;
%     lambda = 16;
%     CapacityCost = [0.2; 0.3];
%     Weights = [1; 0.9];
%     Mu = [0.8;0.6];
% % 
    mu1 = Mu(1,1);
    mu2 = Mu(2,1);

    L = zeros(C1+1, C1+1);
    M1 = zeros(C1+1, C1+1);
    Q = zeros((C1+1)*(C2+1), (C1+1)*(C2+1));
    lowerI = zeros(C1+1, C1+1);
    for i = 1:C1
       L(i, i+1) = lambda;
       M1(i+1,i) = i*mu1;
       lowerI(i+1,i) = 1;
    end
    
    for j = 0:C2
        for i = 1:C1
            Q((C1+1)*j+i, (C1+1)*j+i) = - (lambda+j*mu2);
        end
        Q((C1+1)*j+C1+1, (C1+1)*j+C1+1) = -j*mu2;
        for i = 1:C1
            Q((C1+1)*j+i+1, (C1+1)*j+i+1) = Q((C1+1)*j+i+1, (C1+1)*j+i+1) - i*mu1;
        end
    end
    
    M2 = zeros(C2+1, C2+1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    upperI = zeros(C2+1, C2+1);
    for i = 1:C2
       M2(i+1,i) = i*mu2; 
       upperI(i, i+1) = 1;
    end
    
    D0 = Q + kron(diag([zeros(C2,1); 1]), M1) + kron(eye(C2+1), L) + kron(M2, eye(C1+1));
        
    D1 = kron(upperI, M1);

    D = D0 + D1;
    
    A = D;
    A(:,end) = 1;
    pi = transpose(A)\[zeros((C1+1)*(C2+1)-1,1); 1];
    pi = pi';

    value = pi*(Weights(1,1)*D1)*ones((C1+1)*(C2+1),1)-CapacityCost(1,1)*C1-CapacityCost(2,1)*C2;
end