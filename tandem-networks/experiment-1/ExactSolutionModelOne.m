function [optimalC1,optimalC2,value] = ExactSolutionModelOne(c1max, c2max, Weights, Mu, lambda, CapacityCost)
data = NaN(c1max+1, c2max+1);
%h = waitbar(0,'Finding optimal value using MAMs. This may take more than 10 minutes (and slows down as it progresses).');
for C1 = 0:1:c1max
    for C2 = 0:1:c2max
       data(C1+1,C2+1) = ExactEvaluationModelOne(C1, C2, Weights, Mu, lambda, CapacityCost);
       %waitbar((C2+1+C2max*(C1+1))/((C1max+1)*(C2max+1)),h)
    end
end
%close(h)
% surf(0:1:c2max,0:1:c1max, data),axis('tight'),
% shading('interp'),view([0,90]);
% xlabel('C2'), ylabel('C1')
% c = colorbar;  c.Label.String = 'f(c_1,c_2)';


[value,ind] = max(data(:));
[optimalC1,optimalC2] = ind2sub(size(data),ind);
optimalC1 = optimalC1-1;
optimalC2 = optimalC2-1;
% title(strcat('C_1^*=',num2str(optimalC1),' and C_2^*=',num2str(optimalC2)))
end
