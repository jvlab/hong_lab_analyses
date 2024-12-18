%hlid_connect_eigstats2: compare statistics of eigenvalues of connectivity matrices
% 
%   See also: HLID_MAJAXES, HLID_CSV2CONNECTIVITY_DEMO, HLID_CONNECT_EIGSTATS.
%
if ~exist('fn_connectivity') fn_connectivity='hlid_csv2connectivity_19Sep24.mat'; end
if ~exist('nshuffs') nshuffs=1000; end %number of shuffles to compute
% if ~exist('quantiles_plot') quantiles_plot=[0.99 0.95 0.5 0.05 0.01]; end
if ~exist('if_frozen') if_frozen=1; end %frozen random number gens
if ~exist('neigs_plot') neigs_plot=10; end %number of eigenvalues to plot
if ~exist('nbins') nbins=20; end
%
% get connectivity data
%
fn_connectivity=getinp('file name with connectivity information','s',[],fn_connectivity);
load(fn_connectivity);
connectivity_fields=fieldnames(connectivity);
c_filenames=cell(0);
c_varnames=cell(0);
c_fill_methods=cell(0);
for ic=1:length(connectivity_fields)
    icf=connectivity_fields{ic};
    c_filenames{end+1}=connectivity.(icf).filename;
    c_varnames{end+1}=connectivity.(icf).v_name;
    c_fill_methods{end+1}=connectivity.(icf).fill_method;
end
%
c_filenames=unique(c_filenames);
c_varnames=unique(c_varnames);
c_fill_methods=unique(c_fill_methods);
%
%specify connectivity data
%
% for k=1:length(c_filenames)
%     disp(sprintf('connectivity file %1.0f: %s',k,c_filenames{k}));
% end
for k=1:length(c_varnames)
    disp(sprintf('variable %1.0f: %s',k,c_varnames{k}));
end
c_varname=c_varnames{getinp('choice','d',[1 length(c_varnames)])};
for k=1:length(c_fill_methods)
    disp(sprintf('fill method %1.0f: %s',k,c_fill_methods{k}));
end
c_fill_method=c_fill_methods{getinp('choice','d',[1 length(c_fill_methods)])};
%
ncd=length(c_filenames);
c_desc=cell(ncd,1);
c_data=cell(ncd,1);
for icd=1:ncd
    c_filename=c_filenames{icd};
    c_desc{icd}=sprintf('file: %s var: %s, fill method: %s',c_filename,c_varname,c_fill_method);
    c_string=cat(2,c_filename,'__',c_varname,'__',c_fill_method);
    if isfield(connectivity,c_string)
        c_data{icd}=connectivity.(c_string);
    else
        disp('not found in connectivity file')
    end
    m(:,:,icd)=c_data{icd}.vals_filled;
end
%
m_scale=max(abs(m(:)))*[-1 1];
m_size=size(m,1);
m_pairs=m_size*(m_size-1)/2;
%
m_eivecs=zeros(m_size,neigs_plot,ncd);
for icd=1:ncd
    [eivecs,eivals]=eig(m(:,:,icd));
    [eivals,idx]=sort(diag(eivals),'descend');
    m_eivecs(:,1:neigs_plot,icd)=eivecs(:,idx(1:neigs_plot)); %columns are eigenectors in descending order
end
%compare eigenvectors across the two datasets
hbins=([1:nbins]-0.5)/nbins;
frac=zeros(neigs_plot,neigs_plot);
for icd=1:ncd-1
    for jcd=icd+1:ncd
        dots=m_eivecs(:,:,icd)'*m_eivecs(:,:,jcd);
        if if_frozen
            rng('default');
        end
        dots_shuff=zeros(neigs_plot,neigs_plot,nshuffs);
        for ishuff=1:nshuffs
            dots_shuff(:,:,ishuff)=m_eivecs(:,:,icd)'*m_eivecs(randperm(m_size),:,jcd);
        end
        figure;
        set(gcf,'Position',[100 60 1400 900]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',sprintf('dotprods: ds 1=%s ds2=%s',c_desc{icd},c_desc{jcd}));
        %
        maxy=-Inf;
        for irow=1:neigs_plot
            for icol=1:neigs_plot
                subplot(neigs_plot,neigs_plot,icol+(irow-1)*neigs_plot);
                hist(abs(squeeze(dots_shuff(irow,icol,:))),hbins);
                hold on;
                maxy=max(maxy,max(get(gca,'YLim')));
            end
        end
        for irow=1:neigs_plot
            for icol=1:neigs_plot
                subplot(neigs_plot,neigs_plot,icol+(irow-1)*neigs_plot,'YLim',[0 maxy]);
                plot(repmat(abs(dots(irow,icol)),1,2),[0 maxy],'r','LineWidth',2);
                frac(irow,icol)=sum(abs(dots_shuff(irow,icol,:))>abs(dots(irow,icol)))/nshuffs;
                text(0.4,0.5*maxy,sprintf('p=%5.3f',frac(irow,icol)),'FontSize',7);
                title(sprintf('eig%2.0f.eig%2.0f=%5.3f',irow,icol,dots(irow,icol)));
                set(gca,'XTick',[0 1]);
                set(gca,'YTick',[]);
            end
        end
        %
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,sprintf('connectivity data from %s: set 1: %s, set 2: %s',fn_connectivity,c_desc{icd},c_desc{jcd}),'Interpreter','none');
        axis off;
    end %jcd
end %icd
