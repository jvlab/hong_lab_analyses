function [pc_sel,pc_sel_string,opts_pcasel_new]=hlid_vi_pcaselect(opts_pcasel)
% [pc_sel,pc_sel_string,opts_pcasel_new]=hlid_vi_pcaselect(opts_pcasel)
% selects pcs to include in a volumetric imaging reconstruction
%
% opts_pcasel: options
%   n_pc_max: number of pcs available
%   eiv_squared: (length is n_pc_max) power explained by each pc
%   frats: (length is n_pc_max) f-ratio for excess variance, of across-stim differences
%   fdof: degrees of freedom for frats
%   pcrits: built-in critical p-values, can be omitted
%
%  See also:  HLID_VI_PCAFILT, HLID_VARRATS.
%
opts_pcasel=filldefault(opts_pcasel,'pcrits',[0.001 0.01 0.05]);
%
opts_pcasel_new=opts_pcasel;
%
sel_types_avail={'eigenvalue number or cumulative power','F-ratio for across-stimulus differences'};
sel_types_brief={'power','frat'};
%
if_ok=0;
while (if_ok==0)
    pc_sel=[];
    pc_sel_string=[];
    for is=1:length(sel_types_avail)
        disp(sprintf('%1.0f->%s',is,sel_types_avail{is}));
    end
    sel_list=getinp('1 or more choices (0 will terminate)','d',[0 length(sel_types_avail)]);
    if sel_list(1)==0
        sel_list=[];
        pc_sel_string='none';
    end
    sel_sets=cell(1,length(sel_list));
    sel_meths=cell(1,length(sel_list));
    sel_vals=zeros(1,length(sel_list));
    if length(sel_list)>1
       combine_meth=getinp('1 to combine by AND, 2 by OR','d',[1 2],1);
    else
        combine_meth=2; % works even if sel_list is empty
    end
    switch combine_meth
        case 1
            combine_string='and';
            pc_sel=[1:opts_pcasel.n_pc_max];
        case 2
            combine_string='or';
            pc_sel=[];
    end
    for k=1:length(sel_list)       
        opts_pcasel_new.sel_type{k}=sel_types_avail{sel_list(k)};
        pc_sel_add=[];
        switch sel_types_brief{sel_list(k)}
            case 'power'
                pmeth=getinp('1 by eigenvalue number, 2 by cumulative power','d',[1 2]);
                switch pmeth
                    case 1
                        pc_max=getinp('maximum eigenvalue number','d',[1 opts_pcasel.n_pc_max]);
                        pc_sel_add=sprintf('up to pc %3.0f',pc_max);
                        sel_sets{k}=[1:pc_max];
                        sel_meths{k}='maximum_eigenvalue number';
                        sel_vals(k)=pc_max;
                    case 2
                        vmin_expl=getinp('minimum fraction of variance explained','f',[0 1]);
                        pc_max=max(1,min(find(cumsum(opts_pcasel.eiv_squared)/sum(opts_pcasel.eiv_squared(:))>=vmin_expl)));
                        if vmin_expl>=1
                            pc_max=opts_pcasel.n_pc_max;
                        end
                        pc_sel_add=sprintf('up to pc %3.0f (vmin: %7.5f)',pc_max,vmin_expl);
                        sel_sets{k}=[1:pc_max];
                        sel_meths{k}='minimum fraction variance explained';
                        sel_vals(k)=vmin_expl;   
                end  
            case 'frat'
                fmeth=getinp('1 by min F-ratio, 2 by max p-value','d',[1 2]);
                switch fmeth
                    case 1 %F-ratio
                        frat_min=getinp('minimum F-ratio','f',[0 Inf]);
                        pc_sel_add=sprintf('minimum F-ratio: %8.5f',frat_min);
                        sel_sets{k}=find(opts_pcasel.frats>=frat_min);
                        sel_meths{k}='minimum F-ratio';
                        sel_vals(k)=frat_min;
                    case 2 %p-value
                end
        end
        switch combine_string
            case 'and'
                pc_sel=intersect(pc_sel,sel_sets{k});
            case 'or'
                pc_sel=union(pc_sel,sel_sets{k});
        end
        if k>1
            pc_sel_string=cat(2,pc_sel_string,' ',combine_string,' ',pc_sel_add);
        else
            pc_sel_string=pc_sel_add;
        end
    end
    disp(sprintf('%2.0f pcs chosen (%s), %7.3f of total power',length(pc_sel),pc_sel_string,sum(opts_pcasel.eiv_squared(pc_sel))/sum(opts_pcasel.eiv_squared(:))));
    if_ok=getinp('1 if ok','d',[0 1]);
end
opts_pcasel_new.sel_sets=sel_sets;
opts_pcasel_new.sel_meths=sel_meths;
opts_pcasel_new.sel_vals=sel_vals;
opts_pcasel_new.combine_meth=combine_string;
opts_pcasel_new.sel_types_used=sel_types_avail(sel_list);
return
end
