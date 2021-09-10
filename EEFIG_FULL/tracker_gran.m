%% Tracker gran Function

function [C,m] = tracker_gran(prev_C,prev_m,multiplierN,lambda,xk)
effectiveN = 200;
multiplierN = min([multiplierN effectiveN]);
temp1 = (xk-prev_m)*prev_C*(xk-prev_m)';
temp1 = temp1+(multiplierN-1)/lambda;
mplier = ((multiplierN)/((multiplierN-1)*lambda));
C = prev_C-((prev_C*(xk-prev_m)'*(xk-prev_m)*prev_C)/temp1);
C = mplier*C;
m = lambda*prev_m+(1-lambda)*xk;
end