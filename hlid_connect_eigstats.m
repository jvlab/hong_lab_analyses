%hlid_connect_eigstats: statistics of eigenvalues of connectivity matrices
% 
%   See also: HLID_MAJAXES, HLID_CSV2CONNECTIVITY_DEMO, HLID_CONNECT_EIGSTATS2.
%
if ~exist('fn_connectivity') fn_connectivity='hlid_csv2connectivity_19Sep24.mat'; end
if ~exist('nshuffs') nshuffs=1000; end %number of shuffles to compute
if ~exist('quantiles_plot') quantiles_plot=[0.99 0.95 0.5 0.05 0.01]; end
if ~exist('if_frozen') if_frozen=1; end %frozen random number gens
if ~exist('neigs_sub') neigs_sub=4; end % max number of eigenvectors to subtract
if ~exist('neigs_plot') neigs_plot=10; end %number of eigenvalues to plot
nquantiles=length(quantiles_plot);
lstrings={'orig'};
for iq=1:nquantiles
    lstrings{iq+1}=sprintf('q %.3f',quantiles_plot(iq));
end
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
for k=1:length(c_filenames)
    disp(sprintf('connectivity file %1.0f: %s',k,c_filenames{k}));
end
c_filename=c_filenames{getinp('choice','d',[1 length(c_filenames)])};
for k=1:length(c_varnames)
    disp(sprintf('variable %1.0f: %s',k,c_varnames{k}));
end
c_varname=c_varnames{getinp('choice','d',[1 length(c_varnames)])};
for k=1:length(c_fill_methods)
    disp(sprintf('fill method %1.0f: %s',k,c_fill_methods{k}));
end
c_fill_method=c_fill_methods{getinp('choice','d',[1 length(c_fill_methods)])};
c_desc=sprintf('file: %s var: %s, fill method: %s',c_filename,c_varname,c_fill_method);
c_string=cat(2,c_filename,'__',c_varname,'__',c_fill_method);
if isfield(connectivity,c_string)
    c_data=connectivity.(c_string);
else
    disp('not found in connectivity file')
end
%
m=c_data.vals_filled;
m_scale=max(abs(m(:)))*[-1 1];
m_size=size(m,1);
m_pairs=m_size*(m_size-1)/2;
%
[m_eivecs,m_eivals]=eig(m);
[eivs,idx]=sort(diag(m_eivals),'descend');
m_eivecs=m_eivecs(:,idx); %columns are eigenectors in descending order
%compute eigenvalues of matrices created by shuffling off-diagonal terrms,
%after subtracting some eigenvalues
nrows=2;
ncols=3;
eivals_shuff=zeros(m_size,1+neigs_sub,nshuffs); %d1: eival, d2: number subtracted, d3: shuffle
eivals_shuff_q=zeros(m_size,1+neigs_sub,nquantiles);
m_toshuff=zeros(m_size,m_size,1+neigs_sub);
for isub=0:neigs_sub
    if if_frozen
        rng('default');
    end
    figure;
    set(gcf,'Position',[100 100 1400 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',sprintf('randomized after %2.0f eivecs subtracted',isub));
    subplot(nrows,ncols,1);
    imagesc(m,m_scale);
    axis tight;
    axis square;
    title('orig');
    if (isub==0)
        m_eiv_subtracted=zeros(size(m));
    else
        m_eiv_subtracted=m_eivecs(:,1:isub)*diag(eivs(1:isub))*m_eivecs(:,1:isub)';
    end
    subplot(nrows,ncols,1+ncols);
    imagesc(m_eiv_subtracted,m_scale);
    title(sprintf('eigs up to %2.0f',isub));
    axis tight;
    axis square;
    %
    disp(sprintf(' doing %5.0f shuffles with %2.0f eigenvectors subtracted',nshuffs,isub));
    %
    m_toshuff(:,:,1+isub)=m-m_eiv_subtracted;
    z_toshuff=[]; %harvest upper triangular elements as a single row
    for k=1:m_size-1
        z_toshuff=[z_toshuff,m_toshuff(k,k+1:m_size,1+isub)];
    end
    for ishuff=1:nshuffs
        z=z_toshuff(randperm(m_pairs));
        upper_tri=zeros(m_size);
        zstart=0;
        for k=1:m_size-1 %unpack shuffled elements into a triangle
            upper_tri(k,k+1:m_size)=z(zstart+1:zstart+m_size-k);
            zstart=zstart+m_size-k;
        end
        m_shuffled=m_eiv_subtracted+upper_tri+upper_tri'+diag(diag(m_toshuff(:,:,1+isub)));
        [m_eivecs_shuffled,m_eivals_shuffled]=eig(m_shuffled);
        eivals_shuff(:,1+isub,ishuff)=sort(diag(m_eivals_shuffled),'descend');
        %
        if (ishuff==1)
            isubplot=2;
        elseif (ishuff==nshuffs)
            isubplot=ncols+2;
        else
            isubplot=0;
        end
        if (isubplot>0)
            subplot(nrows,ncols,isubplot)
            imagesc(m_shuffled,m_scale);
            title(sprintf('shuffle %1.0f',ishuff));
            axis tight;
            axis square;
        end
        %also consider shuffles that preserve rowsums and colsums
        %show first and last shuffle in middle column   
    end %ishuff
    %
    eivals_shuff_q(:,1+isub,:)=quantile(eivals_shuff(:,1+isub,:),quantiles_plot,3);
    %
    subplot(1,ncols,ncols);
    plot([1:neigs_plot],eivs(1:neigs_plot),'k.-','LineWidth',2);
    hold on;
    title(sprintf('shuffle after first %1.0f eivecs',isub));
    for iq=1:nquantiles
        if quantiles_plot(iq)>0.5
            qc='r';
        elseif quantiles_plot(iq)<0.5
            qc='b';
        else
            qc='m';
        end
        plot(1:neigs_plot,eivals_shuff_q(1:neigs_plot,1+isub,iq),cat(2,qc,'.:'),'LineWidth',2);
    end
    set(gca,'XTick',[1:neigs_plot]);
    set(gca,'XLim',[0 neigs_plot+0.5]);
    set(gca,'YLim',[0 1.1*eivs(1)]);
    legend(lstrings,'Location','NorthEast');
    hold on;
    %
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,sprintf('connectivity data from %s: %s  (%s)',fn_connectivity,c_string,c_desc),'Interpreter','none');
    axis off;
end
