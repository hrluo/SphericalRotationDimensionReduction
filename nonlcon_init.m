function [c,ceq] = nonlcon_init(x,rad_norm,d_dim)
    %c = sum(abs(x))-rad_norm;
    c = sum(abs(x(1:d_dim) ))-rad_norm;
    ceq = [];
end
%matlab template
%function [c,ceq] = ellipseparabola(x)
%c(1) = (x(1)^2)/9 + (x(2)^2)/4 - 1;
%c(2) = x(1)^2 - x(2) - 1;
%ceq = [];
%end
%
%function [c,ceq] = unitdisk(x)
%c = x(1)^2 + x(2)^2 - 1;
%ceq = [];