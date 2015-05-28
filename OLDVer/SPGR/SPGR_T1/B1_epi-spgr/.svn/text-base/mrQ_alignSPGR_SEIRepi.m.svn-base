function [Res]= mrQ_alignSPGR_SEIRepi(SEIR_T1_1,SEraw1,SEIR_T1_2,sGS2,savepath,saveF,sN,inputAlignFile,OAlignFile)
% Align the SE-IR T1 to SPGR T1
%
% INPUTS:
%   SEIRT1flile : the SE-IR file name with joelle fitting format. 
%                 and a dicom directory of it SEraw1
%                   Magn Reson Med. 2010 Oct;64(4):1057-67. robust 
%                   methodology for in vivo T1 mapping.
%
%   SEIR_T1_2 is the T1 map fitted by multi flip angle on SPGRs
%   SGs is the raw spgr structure we will align all of them
%   savepath path to save
%   Sn name to save
%   inputAlignFile if this thing was done before even partly on this data we
%   can use it (.mat file)
%


%% Check inputs

if (~exist('savepath','var')|| isempty(savepath)),
    savepath=pwd;
end

if (~exist('saveF','var')|| isempty(saveF)),
    saveF=1;
end


%% 

% Load the data
sGS1 = dicomLoadAllSeries(SEraw1);

if exist(OAlignFile,'file')
    load (OAlignFile)
else    
    
    % From here it has to be manual!!!
    answ = questdlg('mrQ_alignSPGR_SEIRepi: This is a user-intensive, manual proccess. Do you want to continue?','mrQ','YES','Get me out of here!','YES');
    if strcmp('YES',answ), answ = 1; else answ = 0; end
    if answ ~=1
        error('User Cancelled: mrQ_alignSPGR_SEIRepi');
    end
    
    edit mrQ_alignSPGR_SEIRepi
    
    if exist(inputAlignFile,'file')
        
        load (inputAlignFile)
        clear Res
        % if you liketo change the mask you got run this
        [f,mn,sd] = defineellipse3d(double(SEIR_T1_1),[],[],mn,sd);
        
        alignvolumedata(SEIR_T1_2,sGS2(1).mmPerVox,SEIR_T1_1,sGS1(1).mmPerVox,ttr);
    else
        
        
        % You need to learn to use knk tool to make a ROI to fit to. See
        % the help defineellipse3d
        [f,mn,sd] = defineellipse3d(double(SEIR_T1_1));
        
        alignvolumedata(SEIR_T1_2,sGS2(1).mmPerVox,SEIR_T1_1,sGS1(1).mmPerVox);
        
        keyboard
        
        % If you got inputAlignFile use this one
        % alignvolumedata(SEIR_T1_2,sGS2(1).mmPerVox,SEIR_T1_1,sGS1(1).mmPerVox,ttr);
        %
        
        % Find  a  good starting point and continue (F5)
        trtemp = alignvolumedata_exporttransformation;
        alignvolumedata_auto(mn,sd,0,[4 4 2]);
        alignvolumedata_auto(mn,sd,0,[2 2 1]);
        alignvolumedata_auto(mn,sd,0,[1 1 1]);
        trRIGID = alignvolumedata_exporttransformation;
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[2 2 1],[],[],[],[],1e-3);
        alignvolumedata_auto(mn,sd,0,[2 2 1]);
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[2 2 1],[],[],[],[],1e-3);
        alignvolumedata_auto(mn,sd,0,[2 2 1]);
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],[],1e-3);
        alignvolumedata_auto(mn,sd,0,[1 1 1]);
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],[],1e-3);
        alignvolumedata_auto(mn,sd,0,[1 1 1]);
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],[],1e-3);
        alignvolumedata_auto(mn,sd,0,[1 1 1]);
        alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],[],1e-3);
        trAFFINE = alignvolumedata_exporttransformation;
    end
    
    keyboard  % if it converge press F5 if not try to do more of the above line until it will
    %play those three untill thre is no improvment!!! if any of those need more
    %then 20 iteration stop it (ctrl C) and go to the next one ...and so on if
    %you can do it and some time run the last four alignvolumedata_auto above.
    %all this is very tricky
    %thing are to slow try the same with [2 2 2]
    % if you are on matlab2008a
    %alignvolumedata_auto(mn,sd,0,[1 1 1]);
    %alignvolumedata_auto(mn,sd,1,[1 1 1]);
    %alignvolumedata_auto(mn,sd,2,[1 1 1]);
    %if you are in matlab9b use this syntaxs. if you.
    %alignvolumedata_auto(mn,sd,0,[1 1 1],[],[],[],[],[],{'Algorithm' ,'lm-line-search'})
    
    % If you are using different version make sure it use the line-search
    % method as this is converge best with it.

    trAFFINE = alignvolumedata_exporttransformation;
    close all;
    ttr=trAFFINE;
end
% [f,mn,sd] = defineellipse3d(double(ll_T1(:,:,:,1)));
% [f,mn,sd] = defineellipse3dAut(double(ll_T1(:,:,:,1)));
% tr=alignvolumedata_exporttransformation;
% alignvolumedata_auto(mn,sd,0,[1 1 1]);

% alignvolumedata_auto(mn,sd,2,[1 1 1]);
% alignvolumedata_auto(mn,sd,3,[1 1 1]);
% alignvolumedata_auto(mn,sd,0,[1 1 1]);
% alignvolumedata_auto(mn,sd,1,[1 1 1]);
% alignvolumedata_auto(mn,sd,0,[1 1 1]);
%ttr=alignvolumedata_exporttransformation;

%% KNK/private

% Change into the private directory of the KNK repo to run the following
% functions: Might want to check for the existence of the KNK repo
% earlier. HACKY
p = pwd; cd(fullfile(mrvDirup(which('alignvolumedata.m')),'private')); %   cd '~avivm/matlab/vistasoft/trunk/kendrick/kendrick/alignvolumedata/private'

ok = reslicevolume(0,ttr,'cubic',3,[],1,0,double(SEIR_T1_2),double(sGS2(1).mmPerVox),size(SEIR_T1_2),double(SEIR_T1_1),double(sGS1(1).mmPerVox),[size(SEIR_T1_1) 1]);

Res{2}.im   = ok;
Res{2}.name = 'align_map';
Res{1}.im   = SEIR_T1_1;
Res{1}.name = 'target_(GS)';

clear ok T1

for i=1:length(sGS2),
    ok=reslicevolume(0,ttr,'cubic',3,[],1,0,double(sGS2(i).imData),double(sGS2(1).mmPerVox),size(SEIR_T1_2),double(SEIR_T1_1),double(sGS1(1).mmPerVox),[size(SEIR_T1_1) 1]);
    Res{i+2}.im=ok;
    Res{i+2}.name=['align_rawFA' num2str(sGS2(i).flipAngle)] ;
    
end

cd (p)

pixdim = sGS1(1).mmPerVox;
imToScanXform = sGS1(1).imToScanXform;


%% Save the results

cd (savepath)
if saveF==1,
    save(sN,'Res' ,'ttr','f','mn','sd', 'pixdim','imToScanXform')
    %save([ sN '/' date 'GSAlign.mat'],'Res' ,'ttr','f','mn','sd', 'pixdim')
end

%plotGSverSPGRB1(Res,1);
