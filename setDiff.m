function [diffVal,diffVec,LHStable,RHStable] = setDiff(LHS,RHS,vars)
% LHS and RHS are created before. This was a choice.
% vars is set1 shift, set1 scale, set2 shift, set2 scale ....
% Determine the number of different sets
stims = [LHS.Properties.RowNames; RHS.Properties.RowNames];

stimID = cell(0);
for stim=stims'
    stimID = [stimID;stim{1}(end-3:end)];
end
stimID = unique(stimID);

% Check vars
%if length(vars) ~= 2*length(stimID)
%    error('wrong number of varaibles');
%end

% Go row by row
LHSmat = LHS{:,:};
RHSmat = RHS{:,:};

for indx = 1:length(stimID)
    % determine variables
    for x=1:length(vars)
        if vars{x}.name == stimID{indx}
            varID = x;
            break;
        end
    end
    moo=contains(LHS.Properties.RowNames,stimID{indx});    
    LHSmat(moo,:) = LHSmat(moo,:)*vars{varID}.scale;
    LHSmat(moo,:) = LHSmat(moo,:)+vars{varID}.shift;
    
    moo=contains(RHS.Properties.RowNames,stimID{indx});
    RHSmat(moo,:) = RHSmat(moo,:)*vars{varID}.scale;
    RHSmat(moo,:) = RHSmat(moo,:)+vars{varID}.shift;

end

RHStable = RHS;
LHStable = LHS;

RHStable{:,:} = RHSmat;
LHStable{:,:} = LHSmat;

diffVec = sqrt(sum((LHSmat-RHSmat).^2,2));
diffVal = sqrt(sum(diffVec.^2));