
% declare your data subject path and subject numbers 
%this line will be different for each user.
% we need to deside what we give with the article:  all data; the result of the loop; an example
% subject,  noting  ect.
%
mapsPath=WH_GetQmrDirs;
sub=[1:19];

%%
% we keep each subject data in stracture WMs{sub}
% we agregate all the subject toghter in a stracture called WM. st is is first point.
st=0;

% loop over subjects
for j=sub
%    get the subject T1 and WF form a mrQ directory organization 
   T1=readFileNifti(fullfile(mapsPath{j}, 'T1_map_lsq.nii.gz'));T1=T1.data;
 WF=readFileNifti(fullfile(mapsPath{j}, 'WF_map.nii.gz'));WF=WF.data;

 
% find each subject segmentation file
FS=readFileNifti(fullfile(fileparts(mapsPath{j}), 'aparc+aseg.nii.gz'));FS=FS.data;
%FS=readFileNifti(fullfile(fileparts(mapsPath{j}), 't1_bet_seg.nii.gz'));FS=FS.data;
%FS=readFileNifti(fullfile(fileparts(mapsPath{j}), 'FS_tissue.nii.gz'));FS=FS.data;

% cheack that the segmentation and the map have the same dimention (some
% time an exsta voxel been added.
SZseg=size(FS); SZdat=size(T1);

if any( SZseg-SZdat~=0);
SZ(1)=min(SZdat(1),SZseg(1))  ;  SZ(2)=min(SZdat(2),SZseg(2))  ;  SZ(3)=min(SZdat(3),SZseg(3))  ;     
    FS=FS(1:SZ(1),1:SZ(2),1:SZ(3));
    T1=T1(1:SZ(1),1:SZ(2),1:SZ(3));
    WF=WF(1:SZ(1),1:SZ(2),1:SZ(3));
end    

% select the WM
%for FSL segmentation
% M= FS==3 & WF<0.88 & WF>0.5; 

% for freee surfare
M=zeros(size(FS));M(FS==2)=1;M(FS==41)=1;M=logical(M);
M= M &WF<0.88 & WF>0.5; 

% calculate where the how many values we have
ed=st+length(find(M==1));
ed1=length(find(M==1));

% agregate the multi subjects data
WM.T1(st+1:ed)=T1(find(M==1));
WM.WF(st+1:ed)=WF(find(M==1));

% save each subject values in WM stracture
WMs{j}.T1(1:ed1)=T1(find(M==1));
WMs{j}.WF(1:ed1)=WF(find(M==1));

st=length(WM.T1);
end



%%
% Bin the R1 values
Xval=[1.05 1.25];
B=linspace(Xval(1),Xval(2),20);


mrvNewGraphWin;
% loop over subjects and plot there R1 vs. 1/WF in each bin R1 bin. 
  for j=sub
      % the mean vaulue for the bin
      MVal=zeros(size(B));
%      the errors of each bin
      SVal =MVal; SEVal=MVal;
      % the values 
      X=1./WMs{j}.T1;
      Y=1./WMs{j}.WF;
      
      for k=1:(length(B)-1) % loop over bins
          % find the location of the R1 in theat bin
          idx=find(X>B(k) & X<B(k+1));
          if ~isempty(idx)
              % calculate the mean STD STE of the bin
              Val=Y(idx);
              MVal(k)=mean(Val);
              SVal(k)=std(Val);
              SEVal(k)=SVal(k)./sqrt(length(Val));
             num(j,k)= length(Val);
          else
              MVal(k)=0;SVal(k)=0;SEVal(k)=0;
          end
      end
      
      % make linear model bitween the R1 and 1/WF bin values 
    idxR=find(MVal>0); y = MVal(idxR); x = B(idxR);
    u = ones(size(idxR));
    V = [x(:), u(:)]; w = inv(V'*V)*V'*y(:)
    
%  plot the sbject mean values and STER and the fitted curve
   plot(x,x.*w(1)+w(2),'k');hold on
      errorbar(B,MVal,SEVal,'k.');
 
 
  end;
  
 ylim([1.35 1.61])
  xlim([1.02 1.3])

 
 
 
%%  repeat for the group data
X=1./WM.T1;
Y=1./WM.WF;

for k=1:(length(B)-1)
          idx=find(X>B(k) & X<B(k+1));
          if ~isempty(idx)
              Val=Y(idx);
              
              MVal(k)=mean(Val);
              SVal(k)=std(Val);
              SEVal(k)=SVal(k)./sqrt(length(Val));
             num(j,k)= length(Val);
          else
              MVal(k)=0;SVal(k)=0;SEVal(k)=0;
          end
      end
 idxR=find(MVal>0); y = MVal(idxR); x = B(idxR);
    u = ones(size(idxR));
    % y = [x u]*w
    V = [x(:), u(:)]; w = inv(V'*V)*V'*y(:)


plot(x,x.*w(1)+w(2),'g')




