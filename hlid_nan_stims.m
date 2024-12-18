%hlid_nan_stims: script to replace one or more stimuli by nans, to set up test files
%
% assumes a coordinate file has been loaded into a structure s
% new structure will be sdel
%
data_fullname=getinp('data file to process','s',[]);
s=load(data_fullname);
nstims=size(s.stim_labels,1);
for istim=1:nstims
    disp(sprintf('%2.0f->stim_label %10s stim_name %15s',istim,s.stim_labels(istim,:),s.stimulus_names(istim,:)));
end
nan_list=getinp('stimuli to replace by nans','d',[1 nstims]);
append_text=cat(2,'_nan',sprintf('-%1.0f',nan_list));
%
dim_prefix='dim';
%
stim_keep=setdiff([1:nstims],nan_list);
snew=s;
snew.dsid=cat(2,s.dsid,append_text);
snew.resps(nan_list(:),:)=NaN;
%
fns=fieldnames(s);
for ifn=1:length(fns)
    if ~isempty(strmatch(dim_prefix,fns{ifn}))
        snew.(fns{ifn})(nan_list(:),:)=NaN;
    end
end
disp('snew');
disp(snew)
data_newname=getinp('new data file to write','s',[],cat(2,data_fullname,append_text));
save(data_newname,'-struct','snew');
disp(sprintf('wrote %s',data_newname));
