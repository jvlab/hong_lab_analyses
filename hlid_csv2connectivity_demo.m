% hlid_csv2connectivity_demo: examine connectivity data
% tablate into a structure "connectivity", along with
% filled-in diagonals, eigenvalues and eigenvectors
%
% connectivity has a field for each file in the library, each variable, and
% each filling-in method
%
% methods for filling in diagonal:
%
%  none: do not fill in
%  zero: fill with zeros
%  eig1: fill in with reconstruction from first eigenvalue and eigenvector
%    (could iterate this to convergence)
%
%  See also: HLID_CSV2CONNECTIVITY.
%
if ~exist('filebases') filebases={'mb','fafb-mb'}; end %library of connectivity data
if ~exist('opts_conn') opts_conn=struct;end
%
fill_methods={'none','zero','eig1'};
plot_vars={'mean','obs';'diff','p';'sd','diff_sd'};
if ~exist('minlog10') minlog10=-3; end %min log for p-values
%
c=struct;
nfiles=length(filebases);
fr=cell(nfiles,1);
ou=cell(nfiles,1);
rd=cell(nfiles,1);
%
for ic=1:length(filebases)
    c_name=strrep(filebases{ic},'-','_');
    [c.(c_name),fr{ic},ou{ic},rd{ic}]=hlid_csv2connectivity(setfield(opts_conn,'file_base',filebases{ic}));
    %
    for iplot=1:size(plot_vars,1);
        figure;
        set(gcf,'Position',[50 100 1400 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,c_name,plot_vars{iplot,1},' and ',plot_vars{iplot,2}));
        %
        for isub=1:size(plot_vars,2)
            subplot(1,size(plot_vars,2),isub);
            v_name=plot_vars{iplot,isub};
            z=c.(c_name).(v_name);
            if strmatch(v_name,'p')
                z=max(log10(z),minlog10);
                z(isnan(z))=minlog10;
                z=-z;
                var_use='-log10(p)';
            else
                var_use=v_name;
            end
            imagesc(z);
            set(gca,'FontSize',7);
            set(gca,'XTick',[1:length(c.(c_name).rois)]);
            set(gca,'XTickLabel',c.(c_name).rois);
            set(gca,'YTick',[1:length(c.(c_name).rois)]);
            set(gca,'YTickLabel',c.(c_name).rois);
            title(var_use,'Interpreter','none');
            axis square;
            colorbar;
        end %subplot
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,sprintf('connectivity data from %s',filebases{ic}),'Interpreter','none');
        axis off;
    end %plot
end
%build a structure of connectivity data, eigenvalues, eigenvectors
connectivity=struct;
for ic=1:length(filebases)
    c_name=strrep(filebases{ic},'-','_');
    v_names=fieldnames(c.(c_name));
    for iv=1:length(v_names)
        v_name=v_names{iv}; 
        if ~strcmp(v_name,'rois')
            z=c.(c_name).(v_name);
            z_orig=z;
            if strmatch(v_name,'p')
                z=max(log10(z),minlog10);
                z(isnan(z))=minlog10;
                z=-z;
                var_use='-log10(p)';
                else
            var_use=v_name;
            end
            disp(' ');
            for im=1:length(fill_methods)
                fill_method=fill_methods{im};
                switch fill_method
                    case 'none'
                        z_filled=z;
                    case 'zero'
                        z_filled=z-diag(diag(z));
                    case 'eig1'
                        z(isnan(z))=0;
                        [v,d]=eigs(z,1); %largest eigenvalue of raw matrix
                        z_recon=v*d*v';
                        z_filled=z-diag(diag(z))+diag(diag(z_recon));
                end
                %
                s_name=cat(2,c_name,'__',v_name,'__',fill_method);
                s=struct;
                s.name=s_name;
                s.filename=c_name;
                s.rois=c.(c_name).rois;
                s.v_name=v_name;
                s.var_use=var_use;
                s.fill_method=fill_method;
                s.vals_orig=z_orig;
                s.vals_transformed=z;
                s.vals_filled=z_filled; %filled values
                s.nancount=sum(isnan(z_filled(:)));
                %eigenvalues and eigenvectors: compute all and then sort
                [v,d]=eig(z_filled);
                d=diag(d);
                [dsort,idx]=sort(d,'descend');
                vsort=v(:,idx);
                s.eigenvalues=dsort;
                s.eigenvectors=vsort;
                %
                connectivity.(s_name)=s;
                disp(sprintf('processed %30s: largest two eigs %10.2f %10.2f  nancount %3.0f',s.name,s.eigenvalues(1:2),s.nancount))
            end %im
        end %not rois
    end %iv
end %ic
disp('suggest saving ''connectivity''');
