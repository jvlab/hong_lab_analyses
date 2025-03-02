function [h,hcol]=tblxinfo_count(tbl_count,lut)
%
% [h,hcol]=tblxinfo_count(tbl_count,lut) computes transinformation and column entropy,
% in bits, from a table of counts.  Saves time by avoiding computation of logs
%
% tbl must be 2-dimensional.
% lut is look-up table of n*log2(n), for n=0 to some large number.  Can be omitted, then lookup table not used.
%
% h is transinformation
% hcol is entropy of column sums
%
% 21Nov20:  allow lut to be omitted
%
% See also TBLXINFO, OPT_TAG_INFO, HISTINFO.
%
if (nargin<=1)
    lut=[];
end
col_count=sum(tbl_count,1);
row_count=sum(tbl_count,2);
ntot=sum(col_count);
if ntot<length(lut)
    hcol=(lut(1+ntot)-sum(lut(1+col_count)))/ntot;
%h=(sum(lut(1+tbl_count(:)))-sum(lut(1+col_count))-sum(lut(1+row_count))+lut(1+ntot))/ntot;
    h=(sum(lut(1+tbl_count(:)))-sum(lut(1+row_count)))/ntot+hcol;
else
    h=tblxinfo(tbl_count/ntot);
    hcol=histinfo(col_count/ntot);
end
return
