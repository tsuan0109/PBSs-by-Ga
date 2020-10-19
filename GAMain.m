function [x, fval, score, bestscore] = GAMain(x0)
%initial structure
%[x, fval, score, bestscore] = GAMain([1.214, 0.929, 0.643, 0.5, 0.5, 0.5, 0.48, 0.48, 0.48, 0.48, 0.48])
score = [];
G = [];
% Initialize
options = gaoptimset(@ga);
% Set Parameters
options = gaoptimset(options, 'CrossoverFraction', 0.7);
options = gaoptimset(options, 'InitialPopulation', x0);
options = gaoptimset(options, 'PopulationSize', 10);
options = gaoptimset(options, 'CreationFcn', @gacreationnonlinearfeasible);
options = gaoptimset(options, 'SelectionFcn',@selectionroulette);
options = gaoptimset(options, 'CrossoverFcn', @crossovertwopoint);
options = gaoptimset(options, 'MutationFcn', @mutationadaptfeasible);
options = gaoptimset(options, 'Display', 'off');
options = gaoptimset(options, 'OutputFcn', @myoutput); 
options = gaoptimset(options, 'TolFun', 1e-1); 
options = gaoptimset(options, 'PlotFcns', {  @gaplotbestf @gaplotbestindiv @gaplotscores });
%Constraint
LB = zeros(1,11);
UB = zeros(1,11);
LB(7:11) = 0.390;
UB(7:11) = 0.570;
LB(1:6) = 0.25;
UB(1:6) = 1.4;
%fitness function
f = @(x) CMT(x);
Conss = @(x) Constraint(x);

[x,fval] = ga(f,11,[],[],[],[],LB,UB,Conss,[],options);

function [state,options,optchanged] = ...
             myoutput(options,state,flag) 
                     B = state.Score;
                     B = B';
          score = [score; B];
          optchanged = false;
end
[bestscore,posi] = min(score');
end
