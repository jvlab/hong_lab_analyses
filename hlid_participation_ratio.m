%hlid_participation_ratio: calculate participation ratio from single-trial data
%
% participation ratio, a measure of effective dimension, is sum(eivs)^2/sum(eivs^2), 
% where eivs are eigenvalues of SVD of response matrix
%
% This works on single-trial datasets, and also computes participation
% ratio based on average across repeats.
%
% Could be adapted to work on trial-averaged datasets too.
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM_TRIAL_READ
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
sub_labels={'',' (mean sub)'}; %subtract mean from responses
nsubs=length(sub_labels);
preproc_labels={'raw','normalized'}; %can replace by a subset to shorten analysis or downsample
npreprocs=length(preproc_labels);
reptavg_labels={'indiv repts','repts avged'};
nreptavgs=length(reptavg_labels);
%
if ~exist('nrepts') nrepts=3; end %number of repeats
%
if ~exist('if_log') if_log=0; end
hlid_rastim_trial_read;
%
% metadata=cell(nfiles,1);
% dsids=cell(nfiles,1);
% resps_mean=cell(nfiles,1); mean responses within stimuli(nstims,nrois_avail)
% trial_ptrs=cell(nfiles,1); array of (nstims,nrepts)
% resps_trial=cell(nfiles,1); trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% trial_sequence=cell(nfiles,1); stimulus sequence, stimuli as strings
% stims_avail=cell(nfiles,1); list of available stimuli in each file, beginning at 1
% rois_avail=cell(nfiles,1); list of roi numbers kept for analysis, beginning at 1
% rois=cell(nfiles,1):  original rois
% nrois_avail(ifile): number of rois available
%
for k=1:nstims
    stimulus_names_display{k}=stimulus_names(k,1:-1+min(find(stimulus_names(k,:)==' ')));
end
partratio=zeros(nfiles,nsubs,npreprocs,nreptavgs); %d1=prep, d2=subtract mean? d3=normalize? d4=avg repts?
eigenvals=cell(nfiles,nsubs,npreprocs,nreptavgs);
%begin revising here, loop over irept_avg and choose resps_mean or resps_trial
for ireptavg=1:nreptavgs
    for ifile=1:nfiles       
        nrois=nrois_avail(ifile);
        switch reptavg_labels{ireptavg}
            case 'indiv repts'
                r=resps_trial{ifile};
                r=reshape(r,[nstims*nrepts,nrois]);
            case 'repts avged'
                r=resps_mean{ifile};
        end %r is now [stim x roi] or [(stim,rept) x roi]
        for isub=1:nsubs
            switch sub_labels{isub}
                case ''
                    ruse=r;
                case ' (mean sub)'
                    rs_xm=mean(r,1,'omitnan'); %global mean
                    ruse=r-repmat(rs_xm,[size(r,1) 1]);
            end
            for ipreproc=1:npreprocs
                switch preproc_labels{ipreproc}
                    case 'raw'
                    case 'normalized'
                        ruse_norm=sqrt(sum(ruse.^2,2));
                        ruse_norm(ruse_norm==0)=1;
                        ruse=ruse./repmat(ruse_norm,[1 nrois]);
                end
                %compute svd
                    nonans_stim=find(~all(isnan(ruse),2)); %exclude stimuli for which all are NaN
                    ruse=ruse(nonans_stim,:);
                    nonans_roi=find(all(~isnan(ruse),1)); %then exclude rois that have any NaN
                    ruse=ruse(:,nonans_roi);
                    %
                    [u,s,v]=svd(ruse); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                    maxrank=min(size(ruse));
                    s=s(1:maxrank,1:maxrank);
                    s=diag(s);
                    eigenvals{ifile,isub,ipreproc,ireptavg}=s;
                    partratio(ifile,isub,ipreproc,ireptavg)=sum(s).^2/sum(s.^2);
            end %ipreproc
        end %isub
        disp(sprintf('processed %14s (%3.0f stims, %4.0f rois, from %3.0f stims, %4.0f rois) for file %s ',...
            reptavg_labels{ireptavg},size(ruse),size(r),dsids{ifile}));
    end %ifile
end %ireptavgs
%
disp(sprintf('statistics of participation ratios for %2.0f files',nfiles));
for ifile=1:nfiles
    disp(sprintf('%40s %40s %s',dsids{ifile},metadata{ifile}.title,metadata{ifile}.imaging_type));
end
disp('                                              mean    min     max     stdv')
for ireptavg=1:nreptavgs
    for isub=1:nsubs
        for ipreproc=1:npreprocs
            label=sprintf('%14s %11s %13s',reptavg_labels{ireptavg},sub_labels{isub},preproc_labels{ipreproc});
            pr=partratio(:,isub,ipreproc,ireptavg);
            disp(sprintf('%30s:  %7.3f %7.3f %7.3f %7.3f',label,mean(pr),min(pr),max(pr),std(pr)))
        end
    end
end
%
%save results
%
results=struct;
results.nrepts=nrepts;
results.nstims=nstims;
results.nfiles=nfiles;
results.nsubs=nsubs;
results.sub_labels=sub_labels;
results.npreprocs=npreprocs;
results.preproc_labels=preproc_labels;
results.nreptavgs=nreptavgs;
results.reptavg_labels=reptavg_labels;
results.metadata=metadata;
results.dsids=dsids;
results.stims_avail=stims_avail; %list of available stimuli in each file, beginning at 1
results.rois_avail=rois_avail;
results.stimulus_names=stimulus_names;
results.stimulus_names_display=stimulus_names_display;
results.partratio=partratio;
results.partratio_dims={'d1: prep, d2: sub mean, d3: normalize, d4: avg repts'};
results.eigenvals=eigenvals;
disp('results structure created.');
