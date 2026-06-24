% analysis may 6

analysis_5_4;

vars = cell(5,1);
vars{1}.name = 'kiTC';
vars{2}.name = 'meTC';
vars{3}.name = 'moSK';
vars{4}.name = 'vaTC';
vars{5}.name = 'wdGB';

% Shift these so that each of the sets has the same original statistics. 


vars{3}.shift = -0.0;
vars{3}.scale = 1.0;

%{
dValMin = 200.0;
varBest = vars;
for iter=1:100000
    for v=1:5
        if(v~=4)
            vars{v}.shift = rand()*6-3;
            vars{v}.scale = rand()*8+1;
        end
    end
    [dVal,dVec,m1,m2] = setDiff(LHS,RHS,vars);
    if(dVal<dValMin)
        dValMin = dVal;
        varBest = vars;
        m1Min = m1;
        m2Min = m2;
    end
        
end
%}

for setindx = 1:5
    if(setindx~=3)
        vars{setindx}.scale = rand()*5;
        vars{setindx}.shift = rand()*10-5;
    end
end

dt=0.001;
dBase = setDiff(LHS,RHS,vars);

grad = zeros(10,1);
varsTmp =  vars;
varsNew = vars;
numIts = 5000;

for iter = 1:numIts
    % Compute the gradient-ish
    for setindx=1:5
        if(setindx ~= 3)
            varsTmp = vars;
            varsTmp{setindx}.shift = vars{setindx}.shift-dt;
            gVal = setDiff(LHS,RHS,varsTmp);
            grad(setindx*2-1) = (dBase-gVal)/dt;
            varsTmp = vars;
            varsTmp{setindx}.scale = vars{setindx}.scale-dt;
            gVal = setDiff(LHS,RHS,varsTmp);
            grad(setindx*2) = (dBase-gVal)/dt;
        end
    end
    
    % Subtract the gradient
    p = [0.001,0.005,0.01,0.1,1.0,10.0];
    cand = zeros(length(p),1);
    for pindx=1:length(p)
        for setindx=1:5
            varsNew{setindx}.shift = vars{setindx}.shift - p(pindx)*grad(setindx*2-1);
            varsNew{setindx}.scale = vars{setindx}.scale - p(pindx)*grad(setindx*2);
        end
        cand(pindx) = setDiff(LHS,RHS,varsNew);
    end
    [val,vix] = min(cand);
    if(val<dBase)
        for setindx=1:5
            varsNew{setindx}.shift = vars{setindx}.shift - p(vix)*grad(setindx*2-1);
            varsNew{setindx}.scale = vars{setindx}.scale - p(vix)*grad(setindx*2);
        end
        dBase = val;
        vars = varsNew;
    else
        break;
    end
    
    
end        
        
    