function h=tblxinfo(tabl)
%
% h=tblxinfo(tabl) gives the transinformation, in bits
% tabl must be 2-dimensional.
%
% sum(tabl) assumed to be 1, and all elements positive.
%
% See also HISTINFO, HISTENT, HISTINFO_NZ, TBLXINFO_COUNT.
%
if any(tabl(:)<=0) %modified 7 Jul 17 to test for any non-positive entries
    h=histinfo(sum(tabl,1))+histinfo(sum(tabl,2))-histinfo(tabl);
else
    h=histinfo_nz(sum(tabl,1))+histinfo_nz(sum(tabl,2))-histinfo_nz(tabl);
end
return

