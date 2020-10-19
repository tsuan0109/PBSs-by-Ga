function [c,ceq] = Constraint(p)
% p = [1.25,1,0.75,0.45,0.45,0.45,0.48,0.48,0.48,0.48,0.48] 
wg = zeros(1,6);
ww = zeros(1,6);
wg(1,1) = p(1)-1.5;
ww(1,1) = p(7)-0.48;
for k=1:4
%wp(1:23) difference between gap  wp(25:48) difference between width of waveguide
    wg(1,k+1) = p(k+1)-p(k);  
    ww(1,k+1) = p(k+7)-p(k+6);
end

wg(1,6) = p(6)-p(5);
ww(1,6) = 0.48-p(11);
wa = zeros(1,6);
wb = zeros(1,6);
for t=1:5
    wa(1,t) = wg(1,t+1)-wg(1,t);
    wb(1,t) = ww(1,t+1)-ww(1,t);
end
    wa(1,6) = wg(1,5)-wg(1,6);
    wb(1,6) = ww(1,5)-ww(1,6);
c(1) = max(abs(wg))-0.44;
c(2) = max(abs(wa))-0.44;

c(3) = max(abs(ww))-0.08; 
c(4) = max(abs(wb))-0.14;
ceq = [];
