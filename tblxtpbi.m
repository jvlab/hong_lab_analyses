function h=tblxtpbi(ctabl,useall)
%
% function h=tblxtpbi(ctabl) gives the Treves-Panzeri bias correctoin
% to be added to the naive transinformation estimate, in bits
%
% ctabl must be 2-dimensional, and contains the counts.
%   sum(ctabl) is the total number of trials
%
% useall=1 to use all bins in the table
% useall=0 (default) to use bins that are in an occupied row and occupied column
% useall=-1 to just use nonzero bins
%
% See also TBLXINFO, HISTTPBI, OPT_SEGMENT_SETUP, OPT_SEGMENT_DEFOPT.
%
if (nargin <=1); useall=0; end
switch useall
	case 1
        h=histtpbi(sum(ctabl,1),1)+histtpbi(sum(ctabl,2),1)-histtpbi(ctabl,1);
	case 0
        pruned=ctabl([find(sum(ctabl,2)>0)],[find(sum(ctabl,1)>0)]);
        h=histtpbi(sum(pruned,1),1)+histtpbi(sum(pruned,2),1)-histtpbi(pruned,1);
   case -1
        h=histtpbi(sum(ctabl,1),0)+histtpbi(sum(ctabl,2),0)-histtpbi(ctabl,0);
end
return
