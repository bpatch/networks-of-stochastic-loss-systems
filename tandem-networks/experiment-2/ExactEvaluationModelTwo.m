function value = ExactEvaluationModelTwo(c1, c2, Weights, Mu, lambda, CapacityCost)
%     c1 = 2;
%     c2 = 3;
%     lambda = 16;
%     CapacityCost = [0.2; 0.3];
%     Weights = [1; 0.9];
%     Mu = [0.8;0.6];
% 
    mu1 = Mu(1,1);
    mu2 = Mu(2,1);

    L = zeros(c1+1, c1+1);
    M1 = zeros(c1+1, c1+1);
    Q = zeros((c1+1)*(c2+1), (c1+1)*(c2+1));
    lowerI = zeros(c1+1, c1+1);
    for i = 1:c1
       L(i, i+1) = lambda;
       M1(i+1,i) = i*mu1;
       lowerI(i+1,i) = 1;
    end
    
    for j = 0:c2
        for i = 1:c1
            Q((c1+1)*j+i, (c1+1)*j+i) = - (lambda+j*mu2);
        end
        Q((c1+1)*j+c1+1, (c1+1)*j+c1+1) = -j*mu2;
        for i = 1:c1
            Q((c1+1)*j+i+1, (c1+1)*j+i+1) = Q((c1+1)*j+i+1, (c1+1)*j+i+1) - i*mu1;
        end
    end
    
    M2 = zeros(c2+1, c2+1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    upperI = zeros(c2+1, c2+1);
    for i = 1:c2
       M2(i+1,i) = i*mu2; 
       upperI(i, i+1) = 1;
    end
    
    LAMBDA = lambda*diag([zeros(c1,1); 1]);
    
    D0 = Q + kron(eye(c2+1), M1) + kron(M2, eye(c1+1)) - kron(diag([ones(c2,1); 0]), LAMBDA);
        
    D1 = kron(eye(c2+1), L);

    D2 = kron(upperI, LAMBDA);
    
    D = D0 + D1 + D2;
    
    A = D;
    A(:,end) = 1;
    pi = transpose(A)\[zeros((c1+1)*(c2+1)-1,1);1];
    pi = pi';

    value = pi*(Weights(1,1)*D1+Weights(2,1)*D2)*ones((c1+1)*(c2+1),1)-CapacityCost(1,1)*c1-CapacityCost(2,1)*c2;
end