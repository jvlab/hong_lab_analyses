function h=histinfo_nz(pvec)
%
% function h=histinfo_nz(pvec) gives the histogram information, in bits
% sum(pvec) assumed to be 1, and all elements positive, but needn't be 1-dimensional.
%
% if some elements may be non-positive, use histinfo.m, which checks and ignores them.
%
% lacks all the bells and whistles of histent.
%
% See also TBLXINFO, HISTENT, HISTINFO.
%
pnz=pvec(:);
h=-pnz'*log(pnz)/log(2);
return

