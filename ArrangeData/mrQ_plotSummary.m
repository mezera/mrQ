function mrQ_plotSummary(mrQ)

%
% This function plots the main output of the mrQ analysis: T1 , T1w, PD (WF), TV,
% SIR Furthermore, it plot important steps in the analysis that might
% indiacte if there is some problem along the way: tissue segmentatio, B1
% and Gain maps, and an overlay of the SEIR and the warpd SPGR to observe
% the registration.

%% open figure
h = figure;set(h, 'Visible', 'off');

%%  go over main maps
% % % % t1
if isfield(mrQ.maps, 'T1path')
    if exist(mrQ.maps.T1path,'file')
        tmp=readFileNifti(mrQ.maps.T1path);
        subplot(5,2,1)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar,
        caxis([0.5,2.5]),         title('T1')
        set(gca,'xTick',[],'yTick',[])
    end
end

% % % % T1w
if isfield(mrQ, 'T1w_files')
    if exist([mrQ.T1w_files,'/T1w.nii.gz'],'file')
        tmp=readFileNifti([mrQ.T1w_files,'/T1w.nii.gz']);
        subplot(5,2,2)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar
        caxis([0,0.2]),         title('T1w')
        set(gca,'xTick',[],'yTick',[])
    end
end

% % % % WF
if isfield(mrQ.maps, 'WFpath')
    if exist(mrQ.maps.WFpath,'file')
        tmp=readFileNifti(mrQ.maps.WFpath);
        subplot(5,2,3)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar
        caxis([0.3,1]),         title('WF')
        set(gca,'xTick',[],'yTick',[])
    end
end

% % % % TV
if isfield(mrQ.maps, 'TVpath')
    if exist(mrQ.maps.TVpath,'file')
        tmp=readFileNifti(mrQ.maps.TVpath);
        subplot(5,2,4)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar
        caxis([0,0.5]),         title('MTV')
        set(gca,'xTick',[],'yTick',[])
    end
end

% % % %  VIP
if isfield(mrQ.maps, 'VIPpath')
    if exist(mrQ.maps.VIPpath,'file')
        tmp=readFileNifti(mrQ.maps.VIPpath);
        subplot(5,2,5)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),       colorbar
        caxis([0,0.5]),        title('VIP')
        set(gca,'xTick',[],'yTick',[])
    end
end
% % % % SIR
if isfield(mrQ.maps, 'SIRpath')
    if exist(mrQ.maps.SIRpath,'file')
        tmp=readFileNifti(mrQ.maps.SIRpath);
        subplot(5,2,6)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar
        caxis([0.2,0.6]),         title('SIR')
        set(gca,'xTick',[],'yTick',[])
    end
end


%% go over bias maps

biadDir=fullfile(mrQ.OutPutNiiDir,'BiasMap');

% % % % B1
if isfield(mrQ, 'B1FileName');
    if exist(mrQ.B1FileName,'file')
        tmp=readFileNifti(mrQ.B1FileName);
        subplot(5,2,7)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),        colorbar
        title('B1')
        set(gca,'xTick',[],'yTick',[])
    end
end

% % % % Gains
Gpath=fullfile(biadDir,'Gains.nii.gz');
if exist(Gpath,'file')
    tmp=readFileNifti(Gpath);
    subplot(5,2,8)
    imagesc(tmp.data(:,:,round(end/2)))
    colormap('gray'),        colorbar
    title('Gains')
    set(gca,'xTick',[],'yTick',[])
end


%% get segmentation file

% % % % Segmentation file
if isfield(mrQ, 'T1w_tissue');
    if exist(mrQ.T1w_tissue,'file')
        tmp=readFileNifti(mrQ.T1w_tissue);
        subplot(5,2,9)
        imagesc(tmp.data(:,:,round(end/2)))
        colormap('gray'),    colorbar,   title('tissue segmentation')
        set(gca,'xTick',[],'yTick',[])
    end
end
%% plot overlay of SEIR and SPGRr

if isfield(mrQ,'SEIR_epi_T1file')  & isfield(mrQ,'Ants_Info')
    if exist(mrQ.SEIR_epi_T1file,'file')  & exist(mrQ.Ants_Info.T1_spgr_epi,'file')
        
        
        
        tmp=readFileNifti(mrQ.SEIR_epi_T1file);
        overIm=readFileNifti(mrQ.Ants_Info.T1_spgr_epi);
        tmp=tmp.data(:,:,round(end/2));overIm=double(overIm.data(:,:,round(end/2)));
        overIm = overIm./max(overIm(:))*256;        tmp = tmp./max(tmp(:))*256;
        
        subplot(5,2,10)
        imagesc((tmp-overIm)./(tmp+overIm))
        colormap('gray'),    colorbar, caxis([-1 1])
        title(sprintf('SEIR-SPGR \n(should have minimal structure)'))
        set(gca,'xTick',[],'yTick',[])
        
        
        % % % not so useful way of overlaying the data
        
        % % % %         Incmap=hot(256);
        % % % %         tmp=tmp.data(:,:,round(end/2));overIm=double(overIm.data(:,:,round(end/2)));
        % % % %         overIm = overIm./max(overIm(:))*256;        tmp = tmp./max(tmp(:))*256;
        % % % %         figure, imagesc(tmp-overIm)  ,colormap('gray'),        colorbar
        % % % %
        % % % %         % Convert to rgb format
        % % % %         h=figure; set(h,'Visible','Off')
        % % % %         colormap(Incmap);  C = colormap;     L = size(C,1);
        % % % %         overIms = round(interp1(linspace(min(overIm(:)),max(overIm(:)),L),1:L,overIm));
        % % % %         overImRgb = reshape(C(overIms,:),[size(overIms) 3]); % Make RGB image from scaled.
        % % % %         close(h)
        % % % %         h=figure; set(h,'Visible','Off')
        % % % %         C = colormap;     L = size(C,1);
        % % % %         tmps = round(interp1(linspace(min(tmp(:)),max(tmp(:)),L),1:L,tmp));
        % % % %         tmpRgb = reshape(C(tmps,:),[size(tmps) 3]); % Make RGB image from scaled.
        % % % %         tmpGray = rgb2gray(tmpRgb);
        % % % %         close(h)
        % % % %
        % % % %
        % % % %         hGray = imshow(tmpGray).   hold on,   hColor = imshow(overImRgb);
        % % % %         set(hColor, 'AlphaData',0.3), title('SPGR over SEIR')  set(gca,'xTick',[],'yTick',[])
    end
    
end


%% make sure it spans the entire "height" of the screen
scrSize=get(0,'ScreenSize');
scrSize(3) = scrSize(3)/2;
set(gcf,'position',scrSize)

imageName1=fullfile(mrQ.OutPutNiiDir,'summary.fig');
imageName2=fullfile(mrQ.OutPutNiiDir,'summary.jpg');

saveas(h,imageName1);
saveas(h,imageName2);
pause(2)
close(h)