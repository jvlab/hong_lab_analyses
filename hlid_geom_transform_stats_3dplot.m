%hlid_geom_transform_stats_3dplot: quick 3-d plot of response space from hlid_gemo transform_stats
%group 1 with solid lines and circles
%group 2 with dotted lines and stars
% frozen random colors used, without affecting the state of the random number gen
% 
% 01Jan26: allow for external creation of title string, for compatibility with hlid_mds_transform_stats
%  See also: HLID_RASTIM_TRIAL_DECODE_3DPLOT, HLID_GEOM_TRANSFORM_STATS, HLID_MDS_TRANSFORM_STATS.
figure;
rng_state=rng;
if (if_frozen~=0) 
    rng('default');
end
if ~exist('title_string')
    title_string=sprintf('iembed %1.0f ->%s, isub %1.0f->%s, ipreproc %1.0f ->%s',iembed,embed_labels{iembed},isub,sub_labels{isub},ipreproc,preproc_labels{ipreproc});
end
%
set(gcf,'Position',[100 100 1400 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',title_string);
colors=rand(nstims,3);
symbs={'o','*'};
%needs work so that one can have multiple points per stimulus if nembed_per_stim>1
if exist('nembed_perstim')
    ieps_lim=nembed_perstim(iembed);
else
    ieps_lim=1;
end
for icol=1:ngps %plot by group
    subplot(1,3,icol)
    for ieps=1:ieps_lim
        for istim=1:nstims
            cplot=reshape(znew_bygp{iembed,idim_ptr}(istim+nstims*(ieps-1),:,gp_list{icol}),[dimlist(idim_ptr),length(gp_list{icol})])'; %prep x dim
            cplot_mean=mean(cplot,1);
            hp=plot3(cplot_mean(1,1),cplot_mean(1,2),cplot_mean(1,3));
            hold on;
            set(hp,'Color',colors(istim,:));
            for iprep=1:size(cplot,1)
                hp=plot3([cplot_mean(1,1) cplot(iprep,1)],[cplot_mean(1,2) cplot(iprep,2)],[cplot_mean(1,3) cplot(iprep,3)],'LineWidth',ieps);
                set(hp,'Color',colors(istim,:));
            end
            hold on;
            title(cat(2,'cons by gp: ',gp_labels{icol},' (',symbs{icol},')'));
            box on;
            grid on;
            axis equal;
            xlabel('d1');
            ylabel('d2');
            zlabel('d3');       
        end
    end
end
%plot groups together, plotting with global consensus
subplot(1,3,3)
for ieps=1:ieps_lim
    for istim=1:nstims
        for icol=1:ngps
            cplot=reshape(znew_glbl{iembed,idim_ptr}(istim+nstims*(ieps-1),:,gp_list{icol}),[dimlist(idim_ptr),length(gp_list{icol})])'; %prep x dim
            cplot_means(icol,:)=mean(cplot,1);
            hp=plot3(cplot_means(icol,1),cplot_means(icol,2),cplot_means(icol,3),symbs{icol});
            hold on;
            set(hp,'Color',colors(istim,:));
        end
        hp=plot3(cplot_means(:,1),cplot_means(:,2),cplot_means(:,3),'LineWidth',ieps);
        set(hp,'Color',colors(istim,:));
        if (ieps==1)
            text(mean(cplot_means(:,1)),mean(cplot_means(:,2)),mean(cplot_means(:,3)),stimulus_names_display{istim});
        end
    end
    title(cat(2,'global consensus (',symbs{1},',',symbs{2},')'));
    box on;
    grid on;
    axis equal;
    xlabel('d1');
    ylabel('d2');
    zlabel('d3');       
end
%
axes('Position',[0.01,0.03,0.01,0.01]); %for text
text(0,0,title_string);
clear title_string
axis off;
rng(rng_state);