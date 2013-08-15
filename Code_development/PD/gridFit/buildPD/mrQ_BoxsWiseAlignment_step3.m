function [M0,Wait,VxUsed,Boxorder,err1]= mrQ_BoxsWiseAlignment_step3(Boxes,ScaleMat,BoxesToUse,boxNm,opt)


if ~isfield(Boxes,'loc')
    Boxes(1).loc=[];
end
if notDefined('boxNm')
    boxNm=1;
end

%% intiate parameters
%brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;
% the box we will work on
DoNow=boxNm;
SZ=size(opt.X);
ToDo=zeros(size(opt.wh));
ToDo(BoxesToUse)= 1;
TodoNext=zeros(size(opt.wh));
doneSt=0;
%out put
M0=zeros(size(BM));
Wait=M0;
count=1;
%%
%loop over boxes
while any (ToDo==1)
    tic
    %for i=1:length(DoNext)
    
    %the over lapping boxes
    wh=find(ScaleMat(:,DoNow)); wh=wh';
    wh=[DoNow wh];
    
    % stack all the over lap voxel in image space
    MMi=nan(length(BM(:)),length(wh));
    Mdc=ScaleMat(wh,DoNow); % this is the scale factor for each box with the box we are working on now
    Mdc(1)=1;
    k=0;
    for ii=wh
        k=k+1;
        if isempty(Boxes(ii).loc)
            XX=Boxes(ii).XX;YY=Boxes(ii).YY;ZZ=Boxes(ii).ZZ;
            mask=zeros(size(BM));
            mask(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2))=1;
            Boxes(ii).loc=find(mask);
        end
        MMi(Boxes(ii).loc(:),k)=Boxes(ii).PD.*Mdc(k);
    end
    
    % find the position where we have at least 6 estimations from differnt
    % boxes
    Im1=sum(~isnan(MMi),2);
    use=find(Im1>6);
    % keep only the position that have multipal estimations
    MMii=nan(length(use),length(wh));
    MMii(:,:)= MMi(use,:);
    clear MMi
    %avrage the estimations
    Im=nanmedian(MMii,2);
    % cheack that the estimation are not wildly different
    ImS=nanstd(MMii,[],2);
    good=(ImS./Im)<0.03;
    Im=Im(good);
    use=use(good);
    % keep track of the number of voxel we add
    VxUsed(DoNow,1)=length(use);
    
    % add the new voxel to the great M0 image
    if doneSt==0 % cheack if this is the first time we add voxels
        doneSt=1;
        M0(use)=Im;
        Wait(use)=1;
        VxUsed(DoNow,2)=VxUsed(DoNow,1);
        err1(DoNow)=0;
    else
        whV=find(M0);
        % we need to find the overlap btween this voxels and the great M0.
        [tf, loc]=ismember(use,whV); % the overlap locations
        loc=loc(loc>0);
        if ~isempty(loc)
            % we need to calculate the offset btween the added voxel and the great M0.
            dc=median(M0(whV(loc))./Im(tf));
            VxUsed(DoNow,2)=length(find(tf==0));
            
            err1(DoNow)=median(abs (M0(whV(loc))-Im(tf).*dc)./M0(whV(loc)));
            %  lets be sure that this is not totly off
            if (err1(DoNow)<0.01 && dc>0)
                % we adding this voxel as wited sum. in the location to the location we
                % allready estimate.
                M0(use)=(M0(use).*Wait(use)+Im.*dc)./(Wait(use)+1); %waited sum
                Wait(use)=Wait(use)+1;
            end
        end
    end
    %book keeping
    Boxorder(DoNow)=count;
    count=count+1
    ToDo(DoNow)=0;
    clear Im1 MMii use Im good ImS Mdc tf loc dc whV
    
    %% get ready for the next round
    % what box we will do next?
    %we can try any of the box that was overlap this box and was not done
    %yet
    TodoNext(wh(2:end))=1;
    TodoNext(ToDo==0)=0;
    % clear all the vriable in the while loop
    
    %  find the box that  maxized overlap to the great M0 defined voxels
    wh=find(TodoNext);
    whV=find(M0);
    for ii=1:length(wh)
        tf =ismember(Boxes(ii).loc(:),whV); % the overlap locations
        overlap(ii)=length(find(tf));
    end
    % if none of the candidate box overlap the great M0 we might be done.
    % but we can try all the other box that are left to be sure
    if notDefined('overlap')
        % all the box that was not tryied yet
        toTry=find( ToDo==1);
        if isempty(toTry) % in this case we are done
            ToDo(:)=0;
        else
            DoNow=toTry(1); % if not we go on and try 
        end
    else
        % next we will work with the box that overlap the great M0 the
        % most.
        mostOverLap=find(overlap==max(overlap));
        DoNow=wh(mostOverLap(1));
        clear wh whV tf  overlap
        toc
    end
    
end
A=1;



