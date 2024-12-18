%hlid_remove_stims: script to remove one or more stimuli, to set up test files
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
remove_list=getinp('stimuli to remove','d',[1 nstims]);
append_text=cat(2,'_deleted',sprintf('-%1.0f',remove_list));
%
dim_prefix='dim';
%
stim_keep=setdiff([1:nstims],remove_list);
snew=struct;
snew.dsid=cat(2,s.dsid,append_text);
snew.metadata=s.metadata;
snew.coord_opts=s.coord_opts;
snew.resps=s.resps(stim_keep(:),:);
snew.stimulus_names=s.stimulus_names(stim_keep(:),:);
snew.stim_labels=s.stim_labels(stim_keep(:),:);
%
fns=fieldnames(s);
for ifn=1:length(fns)
    if ~isempty(strmatch(dim_prefix,fns{ifn}))
        x=s.(fns{ifn});
        snew.(fns{ifn})=x(stim_keep(:),:);
    end
end
disp('snew');
disp(snew)
data_newname=getinp('new data file to write','s',[],cat(2,data_fullname,append_text));
save(data_newname,'-struct','snew');
disp(sprintf('wrote %s',data_newname));
