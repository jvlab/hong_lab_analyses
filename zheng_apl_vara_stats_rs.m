%zheng_apl_vara_stats_rs:  statistics of comparison of subsets, TNT and control files 
% Modular version of zheng_apl_vara_stats
%run after zheng_apl_read_align_rs.m
%
%   See also: RS_VARA_COORDSETS.
% 
%customize to change the quantiles shown
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
%
if ~isfield(preptypes,data_use)
    error('These files have only one kind of prep and is not suitable for this analysis.');
end
if_savemat=getinp('1 to save key variables in mat-file','d',[0 1]);
if_savefig=getinp('1 to save figure(s)','d',[0 1]);
file_name_base=cat(2,'zheng_apl_vara_stats_',data_use,sm_string);
if (if_savemat | if_savefig)
    if_ok=0;
    while (if_ok==0)
        file_name_base=getinp('file name base','s',[],file_name_base);
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
aux=struct;
opts_vara=struct();
opts_vara.allow_reflection=1;
opts_vara.allow_offset=1;
opts_vara.allow_scale=0;
opts_vara.if_normscale=1;
opts_vara.max_niters=1000;
opts_vara.pcon_init_method=0;
opts_vara.if_stats=1;
opts_vara.if_frozen=1;
opts_vara.if_log=1;
%
opts_vara.if_exhaust=0;
if if_debug
    opts_vara.nshuffs=10;
else
    opts_vara.nshuffs=500;
end
opts_vara.shuff_quantiles=shuff_quantiles;
pcon_dim_max=nstims_all-dim_reduce;
max_dim_all=pcon_dim_max;
opts_vara.dim_list_in=[1:pcon_dim_max];
%
aux.opts_vara=opts_vara;
groupings=struct();
groupings.gps=preptypes.(data_use);
[vara_stats,aux_out]=rs_vara_coordsets(data_aligned,groupings,aux)
%
if (if_savefig)
    savefig(gcf,file_name_fig);
    disp(sprintf('figure saved as %s',file_name_fig));
end
%
%save mat file
%
if (if_savemat)
    file_name_mat=file_name_base;
    save(file_name_mat,'data_disp*','opts*','-v7.3');
    disp(sprintf('variables saved as %s',file_name_mat));
end
