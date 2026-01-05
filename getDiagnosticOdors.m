function Sout = getDiagnosticOdors(ref_stimuli,S)

Sout = S;

% remove the last four characters from each of the reference stimuli.
for stimindx = 1:length(ref_stimuli)
    rss(stimindx) = strip(string(ref_stimuli(stimindx,:)));
    rss(stimindx) = extractBefore(rss(stimindx),strlength(rss(stimindx))-4+1);
end

rss = unique(rss);

% Go through the sets and determine the matching diagnostic odors for each.

kept_stims = {};


odorList = cell(0);
setStim = S.Properties.RowNames;
numStimsToCheck = length(setStim);
numRefStims = length(rss);
for refStim = 1:numRefStims
    for checkStims=1:numStimsToCheck
        %rss(refStim)
        %setStim{checkStims}
        if(contains(setStim{checkStims},rss(refStim)))
            odorList = [odorList, setStim{checkStims}];
        end
    end
end
    %numFiles = length(S{setindx});
    %for fileindx = 1:numFiles
    %    Sout{setindx}{fileindx}=S{setindx}{fileindx}(odorList,:);
    %end
odorList
%S(odorList(:),:)

end
