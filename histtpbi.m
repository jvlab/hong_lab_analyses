function h=histtpbi(cvec,useall)
% h=histtpbi(cvec,useall) gives the Treves-Panzeri bias that should be added
% to the naive plugin estimate, in bits
%
% cvec is the vector of counts
% useall=0 (default) to just count the nonzero bins in pvec; 1 to use all bins
%
% See also TBLXINFO, TBLXTPBI, HISTENT, HISTINFO, HISTJABI.
%
if (nargin <=1); useall=0; end
counts=reshape(cvec,1,prod(size(cvec)));
if (sum(counts)==0)
   h=0;
   return
end
if (useall)
   bins=length(counts);
else
	bins=sum(counts>0);
end
h=(bins-1)/(2*sum(counts)*log(2));
return

