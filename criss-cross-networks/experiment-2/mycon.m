function [c,ceq] = mycon(x)
    c = -ones(length(x),1);     % Compute nonlinear inequalities at x.
    ceq = zeros(length(x));   % Compute nonlinear equalities at x.
end