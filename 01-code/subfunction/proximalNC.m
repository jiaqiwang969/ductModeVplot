function [x] = proximalNC(v, lambda, zeta)

% This function computes the proximal operator of a weakly convex J penalty.
% It solves the optimization
%         argmin_{x} lambda * J_{zeta}(x)+(1/2)||x-v||_2^2
% where J_{zeta}(x)=\sum_{i=1}^N F_{zeta}(x_i) and F_{zeta}(x_i) is defined as
%         F_{zeta}(x_i) = |x_i|-zeta*x*x      |x_i| <= 1/2/zeta
%                       = 1/4/zeta            |x_i| > 1/2/zeta
%
% Written by Laming Chen (chen-lm06@mails.tsinghua.edu.cn).

if 2*lambda*zeta >= 1
    x = (2*zeta*abs(v).^2>lambda).*v;
else
    x = (abs(v)>lambda).*(2*zeta*abs(v)<=1).*((v-lambda*sign(v))/(1-2*lambda*zeta))...
        +(2*zeta*abs(v)>1).*v;
end