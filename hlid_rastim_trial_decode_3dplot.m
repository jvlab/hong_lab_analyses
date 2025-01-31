%hlid_rastim_trial_decode_3dplot: quick 3-d plot of response space
%in-sample with solid lines and circles
%out-of-sample with dotted lines and stars
% frozen random colors used, without affecting the state of the random number gen
% 
%  See also: HLID_RASTIM_TRIAL_DECODE.
figure;
rng_state=rng;
if (if_frozen~=0) 
    rng('default');
end
set(gcf,'Position',[100 100 1400 800]);
for istim=1:nstims
    colors=rand(1,3);
    coords_in_plot=coords_insample_align{istim};
    coords_out_plot=coords_outsample_align{istim};
    hp=plot3(coords_in_plot(:,1),coords_in_plot(:,2),coords_in_plot(:,3),'o');
    set(hp,'Color',colors);
    set(hp,'MarkerSize',10);
    hold on;
    coords_in_plot_mean=mean(coords_in_plot,1,'omitnan');
    text( coords_in_plot_mean(1), coords_in_plot_mean(2), coords_in_plot_mean(3),stimulus_names_display{istim});
    if ~isempty(coords_out_plot)
        hp=plot3(coords_out_plot(:,1),coords_out_plot(:,2),coords_out_plot(:,3),'kh');
        set(hp,'Color',colors);
        set(hp,'MarkerSize',10);
    end
    %
    for ip=1:size(coords_in_plot,1)
        hp=plot3([coords_in_plot(ip,1) coords_in_plot_mean(1)],[coords_in_plot(ip,2) coords_in_plot_mean(2)],[coords_in_plot(ip,3) coords_in_plot_mean(3)],'k');
        set(hp,'Color',colors);
        set(hp,'LineWidth',1);
    end
    for ip=1:size(coords_out_plot,1)
        hp=plot3([coords_out_plot(ip,1) coords_in_plot_mean(1)],[coords_out_plot(ip,2) coords_in_plot_mean(2)],[coords_out_plot(ip,3) coords_in_plot_mean(3)],'k:');
        set(hp,'Color',colors);
        set(hp,'LineWidth',2);          
    end
end
box on;
grid on;
axis equal;
xlabel('d1');
ylabel('d2');
zlabel('d3');
%
axes('Position',[0.01,0.08,0.01,0.01]); %for text
text(0,0,sprintf('subsamp %2.0f of %2.0f',isubsamp,nsubsamps_use),'Interpreter','none');
axis off;
%
axes('Position',[0.01,0.05,0.01,0.01]); %for text
text(0,0,sprintf('%s %s ixv_make %3.0f ifold %3.0f',sub_labels{isub},preproc_labels{isub},ixv_make,ifold),'Interpreter','none');
axis off;
%
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,xv_label,'Interpreter','none');
axis off;
rng(rng_state);