function [c,files_read,opts_used,rawdata]=hlid_csv2connectivity(opts)
% [c,files_read,opts_used,rawdata]=hlid_csv2connectivity(opts) reads a set of connectivity files
%
% opts.file_base: base file name, if empty, 'mb' (also could be 'fafb-mb')
% opts.file path: path, defaults to 'C:\Users\jdvicto\Dropbox\From_HongLab\HongLabOrig_for_jdv\connectivity' if empty
% opts.file_list: defaults to {'-matrix','-nm-stats','-result'} ('-nm-counts' omitted)
%
% c: structure of connectivity data
%    c.rois: names of rois, as a cell array
%    other fields are the fields of the csv files:
%      -nm-stats yields mean, sd;
%      -result yields obs, diff, diff-sd and p
% files_read: list of files successfully read
% opts_used: options used
% rawdata: structure of fields
% 
% full file name is file_path\file_base-file_list{k}
%
%readme from Betty Hong:
% ***PN-KC pairwise convergence analysis for FAFB and hemibrain datasets***
% 
% - mb/*
%     Hemibrain MB data.
% - fafb-mb/*
%     FAFB MB data.
% - x/x-matrix.csv
%     Observed connectivity matrix.
% - x/x-nm-counts.csv
%     Co-draw counts for each pair of glomeruli in each of 10,000 shuffled null models.
% - x/x-nm-stats.csv
%     Co-draw count statistics for the null models.
% - x/x-result.csv
%     Pairwise convergence z-scores and p-values for the observed matrix.
%  See also:  HLID_MAJAXES, HLID_CSV2CONNECTIVITY_DEMO, FILLDEFAULT.
%
if nargin<1
    opts=struct;
end
opts=filldefault(opts,'file_base','mb');
opts=filldefault(opts,'file_path','C:\Users\jdvicto\Dropbox\From_HongLab\HongLabOrig_for_jdv\connectivity');
opts=filldefault(opts,'file_list', {'-matrix','-nm-stats','-result'});
%
nfiles=length(opts.file_list);
rawdata=struct;
files_read=cell(nfiles,1);
for ifile=1:nfiles
    fullname=cat(2,opts.file_path,filesep,opts.file_base,opts.file_list{ifile},'.csv');
    if ~exist(fullname,'file')
        disp(sprintf('cannot find %s',strrep(fullname,'\','/')));
    else
        fieldname=strrep(opts.file_list{ifile},'-','_');
        if fieldname(1)=='_'
            fieldname=fieldname(2:end);
        end
        rawdata.(fieldname)=readcell(fullname);
        files_read{ifile}=fullname;
        disp(sprintf('read %s',strrep(fullname,'\','/')));
    end
end
%accumulate all ROI names
rois=cell(0);
table_fields={'nm_stats','result'};
for it=1:length(table_fields)
    tf=table_fields{it};
    if isfield(rawdata,tf)
        rois_new=unique([rawdata.(tf)(2:end,1);rawdata.(tf)(2:end,2)]); %first two cols are an roi pair
        rois=unique([rois_new;rois_new]);
        disp(sprintf('%3.0f unique rois found in %12s, %3.0f rois so far',length(rois_new),tf,length(rois)));
    end
end
nrois=length(rois);
%
c=struct;
c.rois=rois;
for it=1:length(table_fields)
    tf=table_fields{it};
    if isfield(rawdata,tf)
        inds=zeros(size(rawdata.(tf),1)-1,2);
        match_none=0;
        match_mult=0;
        for ientry=1:size(rawdata.(tf),1)-1
            for ij=1:2
                imatch=strmatch(rawdata.(tf){ientry+1,ij},rois,'exact');
                if length(imatch)==0
                    match_none=match_none+1;
                elseif length(imatch)>1
                    match_mult=match_mult+1;
                else
                    inds(ientry,ij)=imatch;
                end
            end
        end
        disp(sprintf(' number of rois in %10s without match: %2.0f, multiple matches: %2.0f',tf,match_none,match_mult));
        vars=rawdata.(tf)(1,3:end);
        for iv=1:length(vars)
            vn=vars{iv};
            c.(vn)=NaN(nrois,nrois);
            for ientry=1:size(rawdata.(tf),1)-1
                if all(inds(ientry,:))>0
                    z=rawdata.(tf){ientry+1,2+iv};
                    if isnumeric(z) & length(z)==1
                        c.(vn)(inds(ientry,1),inds(ientry,2))=z;
                    end
                end
            end %ientry
            n_missing=sum(isnan(c.(vn)(:)));
            disp(sprintf('processed %10s in %10s: %5.0f entries filled in, %5.0f entries not filled in',vn,tf,prod(size(c.(vn)))-n_missing,n_missing));
        end
    end
end
%
opts_used=opts;
