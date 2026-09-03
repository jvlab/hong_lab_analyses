function meths=hlid_meths_define(opts)
% meths=hlid_meths_define(opts) defines methods for reduction of response ectors to coordinates 
%
%  See also: HLID_RASTIM_MDS_COORDS_MAKE, HLID_RASTIM_MDS_COORDS_DEMO, HLID_VI_PCAFIT_COORDS_AUTO.
%
if nargin<1
    opts=struct;
end
meths=cell(0);
meths{1}.name_full='Euclidean distance via SVD';
meths{1}.name_short='Euc dist SVD';
meths{1}.xform='none';
meths{1}.dimred='svd';
meths{1}.name_file='euc_svd';
%
meths{2}.name_full='Euclidean distance via MDS';
meths{2}.name_short='Euc dist MDS';
meths{2}.xform='none';
meths{2}.name_file='euc_mds';
%
meths{3}.name_full='cosine similarity'; %1-normalized dot product
meths{3}.name_short='cos sim';
meths{3}.xform='1-dp';
meths{3}.name_file='cos_sim';
%
meths{4}.name_full='cosine similarity as angle';
meths{4}.name_short='cos ang';
meths{4}.xform='acos(dp)';
meths{4}.name_file='cos_ang';
%
meths{5}.name_full='cosine similarity as chord';
meths{5}.name_short='cos chord';
meths{5}.xform='sqrt(2)*sqrt(1-dp)';
meths{5}.name_file='cos_chord';
%
meths{6}.name_full='Pearson similarity'; %1-normalized centered dot product
meths{6}.name_short='Pearson sim';
meths{6}.xform='1-dp';
meths{6}.name_file='pears_sim';
%
meths{7}.name_full='Pearson similarity as angle';
meths{7}.name_short='Pearson ang';
meths{7}.xform='acos(dp)';
meths{7}.name_file='pears_ang';
%
meths{8}.name_full='Pearson similarity as chord';
meths{8}.name_short='Pearson chord';
meths{8}.xform='sqrt(2)*sqrt(1-dp)';
meths{8}.name_file='pears_chord';
%
return
end
