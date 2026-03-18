function [stims,opts_used]=hlid_vi_stimnames(opts)
% [stims,opts_used]==hlid_vi_stimnames(opts) provides stimulus names and display orders for
% volumetric imaging data from George Barnum, Hong Lab
%
% opts fields:
%  trunc_length: length to truncate names_short
%
% stims.names: stimulus names, as a cell array\
% stims.names_short: short names, abbevs shosen to match other datasets
% stims.display_order: display order
%
% opts_used: options used
%
%   See also:  HLID_VI_EXPLORE,HLID_VI_READ, HLID_SETUP_FUNC.
%
if nargin<1
    opts=struct;
end
opts=filldefault(opts,'trunc_length',6);
opts_used=opts;

%
% From G. Barnmum: Stimulus ordering in the hdf5 and the stimulus identity csv are the same, and the provided sort index csv is a list of indexes into them
% (e.g. the first entry of gbarnum_sort_index_1_indexed.csv is 24, so the 24th entry in the data file and the stimulus identity (CO_2) should be displayed 1st ).
%
stims=struct;
%from  gbarnum_odor_names.csv:
stims.names={...
    'p-cresol',...
    'geosmin',...
    'isoamyl acetate',...
    'banana',...
    'menthone',...
    'mint',...
    '3-methylthio-1-propanol',...
    '1,4 diaminobutane',...
    '2,5-dimethylpyrazine',...
    'geranyl acetate',...
    'trans-2-hexenal',...
    'methyl salicylate',...
    'pfo',...
    '4-methylcyclohexanol',...
    '3-octanol',...
    '1-hexanol',...
    'methyl acetate',...
    '2-heptanone',...
    'acetic acid',...
    'ACV 10%',...
    'ethanol 15%',...
    'wine',...
    'water',...
    'CO2',...
    };
stims.names_short=stims.names;
stims.names_short=strrep(stims.names_short,'1,4 diaminobutane','1,4d');
stims.names_short=strrep(stims.names_short,'1-hexanol','1-6ol');
stims.names_short=strrep(stims.names_short,'2,5-dimethylpyrazine','2,5d');
stims.names_short=strrep(stims.names_short,'2-heptanone','2h');
stims.names_short=strrep(stims.names_short,'3-methylthio-1-propanol','3mt1p');
stims.names_short=strrep(stims.names_short,'3-octanol','3-8ol');
stims.names_short=strrep(stims.names_short,'4-methylcyclohexanol','4-mch');
stims.names_short=strrep(stims.names_short,'ACV 10%','ACV');
stims.names_short=strrep(stims.names_short,'CO2','CO2');
stims.names_short=strrep(stims.names_short,'acetic acid','aa');
stims.names_short=strrep(stims.names_short,'banana','ban');
stims.names_short=strrep(stims.names_short,'ethanol 15%','EtOH');
stims.names_short=strrep(stims.names_short,'geosmin','geos');
stims.names_short=strrep(stims.names_short,'geranyl acetate','ga');
stims.names_short=strrep(stims.names_short,'isoamyl acetate','IaA');
stims.names_short=strrep(stims.names_short,'menthone','menth');
stims.names_short=strrep(stims.names_short,'methyl acetate','ma');
stims.names_short=strrep(stims.names_short,'methyl salicylate','ms');
stims.names_short=strrep(stims.names_short,'mint','mint');
stims.names_short=strrep(stims.names_short,'p-cresol','p-cre');
stims.names_short=strrep(stims.names_short,'pfo','pfo');
stims.names_short=strrep(stims.names_short,'trans-2-hexenal','t2h'); 
stims.names_short=strrep(stims.names_short,'water','H2O'); 
stims.names_short=strrep(stims.names_short,'wine','wine'); 
for k=1:length(stims.names_short)
    if length(stims.names_short{k})>opts.trunc_length
        stims.names_short{k}=stims.names_short{k}(1:opts.trunc_length);
    end
end
%
%from gbarnum_sort_index_1_indexed.csv:
stims.display_order=[24	1	5	3	8	6	15	16	13	12	18	9	23	2	17	11	10	7	19	14	21	4	22	20];
opts=filldefault(opts,'if_log',1);
opts_used=opts;
%
return
