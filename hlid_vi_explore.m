%initial exploration of KC volumetric imaging, data from George Barnum, Hong Lab
%
%  Need to add: plot as function of time
% 
%   See also:  HLID_VI_READ, HLID_VI_STIMNAMES.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('nstims') nstims=24; end
if ~exist('nrepts') nrepts=5; end
%
resp_measures={'deltaF/F','z'};
for k=1:length(resp_measures)
    disp(sprintf('%1.0f->response measure %s',k,resp_measures{k}));
end
resp_measure=resp_measures{getinp('choice','d',[1 length(resp_measures)],1)};
%
if ~exist('if_reorder') if_reorder=0; end
%
rept_list=getinp('repeat list','d',[1 nrepts],1:nrepts);
stim_list=getinp('stimulus list','d',[1 nstims],1:nstims);
opts_read.data_path=data_path;
opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
[s,opts_read_used]=hlid_vi_read(opts_read);
%
disp(s)
%
stims=hlid_vi_stimnames;
%
%olot mean response (delta-F/F), averaged over repeats
%
resp_maxlength=size(s.responses,2);
deltaF=s.responses-repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
switch resp_measure
    case 'deltaF/F'
        v=deltaF./repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
    case 'z'
        v=deltaF./repmat(reshape(s.baseline_stdvs,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
end
v_range=[min(v(:)) max(v(:))];
v_perstim=reshape(mean(mean(v,2,'omitnan'),3,'omitnan'),[s.n_pixels_kept,s.n_stims_kept]); %average across time and repeat
%
%make a heatmap of each slice, for each stimulus
%
v_perstim_range=[min(v_perstim(:)) max(v_perstim(:))];
xyz_range=double([min(s.xyz_kept);max(s.xyz_kept)]);
x_len=xyz_range(2,1)-xyz_range(1,1)+1;
y_len=xyz_range(2,2)-xyz_range(1,2)+1;
[nrows,ncols]=nicesubp(s.n_stims_kept,0.7);
xyz_rel=double(s.xyz_kept)-repmat(xyz_range(1,:),s.n_pixels_kept,1)+1; %relative index of all kept pixels
%
display_ptr_order=[1:s.n_stims_kept];
if if_reorder
    display_sort=zeros(1,s.n_stims_kept);
    for k=1:s.n_stims_kept
        display_sort(k)=find(stims.display_order==s.opts_read.stim_list(k));
    end
    [dsort,display_ptr_order]=sort(display_sort);
end
for plane_ptr=1:s.n_planes_with_data_kept
    plane=s.plane_list_kept(plane_ptr);
    pxls_inplane=find(s.xyz_kept(:,3)==plane);
    tstring=sprintf('plane %2.0f: %s, %s',plane,resp_measure,strrep(s.opts_read.data_file,'.hdf5',''));
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring);
    for stim_ptr_seq=1:s.n_stims_kept
        stim_ptr=display_ptr_order(stim_ptr_seq);
        stim_no=opts_read.stim_list(stim_ptr);
        stim_name=stims.names_short{stim_no};
        heatmap=NaN(x_len,y_len);
        for ipxl=pxls_inplane'
            heatmap(xyz_rel(ipxl,1),xyz_rel(ipxl,2))=v_perstim(ipxl,stim_ptr);
        end
        subplot(nrows,ncols,stim_ptr_seq);
        imagesc(heatmap',v_perstim_range); %matlab plots the transpose
        axis equal;
        set(gca,'XTick',[1 x_len]);
        set(gca,'XTickLabels',xyz_range(:,1));
        set(gca,'YTick',[1 y_len]);
        set(gca,'YTickLabels',xyz_range(:,2));
        axis tight;
        colorbar;
        title(sprintf('s %2.0f: %s',stim_no,stim_name));
    end
    axes('Position',[0.01,0.01,0.01,0.01]);
    text(0,0,tstring,'Interpreter','none');
    axis off
end
