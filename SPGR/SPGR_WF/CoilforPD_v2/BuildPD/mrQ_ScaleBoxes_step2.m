function [Boxes, scaleFactor]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt,errTres,BMfile)
%
% [Boxes, scaleFactor]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt,errTres,BMfile)
%
% This is step 2 of 6 for building the WF map.
%     Step 2: Find the scalar between each box's PD      
%
% ~INPUTS~ 
%        Boxes:
%   BoxesToUse:
%          opt:
%      errTres: Default is 0.01.
%       BMfile:
%
% ~OUTPUTS~
%        Boxes:
%  scaleFactor:
%
% See also: mrQ_buildPD_ver2
%
% AM (C) Stanford University, VISTA

if notDefined('errTres')
        % default for any other case (include T1 reg) 1% error
    errTres=0.01;

%     % the max median present error between two boxes that we still combine.
%     if isfield(opt,'T1reg');
%         % in case we don't use T1 regularization themedian present error defult
%         % treshold is higher. (5%)
%         if opt.T1reg==0
%             errTres=0.05;
%         end
%         
%     end

end

kk=0;

%LinScaleMat=zeros(32*length(Boxes), length(Boxes));
%LinScaleMat=zeros(length(Boxes), length(Boxes));
scaleFactor=zeros(length(Boxes), length(Boxes));

donemask=ones(length(Boxes), length(Boxes));
donemask(BoxesToUse,BoxesToUse)=0;
%Ref=randperm(BoxesToUse);
if notDefined('BMfile')
    BMfile=opt.BMfile;
end

BM=readFileNifti(BMfile);
BM=BM.data;

Boxes(1).loc=[];
for ii=BoxesToUse %loop over boxes
   % tic
   
    donemask(ii,ii)=1;
    % What are the candidate boxes to overlap with the ii-th box?
    wh=find(donemask(ii,:)==0);
    % the ii box's location
    XX=Boxes(ii).XX;YY=Boxes(ii).YY;ZZ=Boxes(ii).ZZ;
    
    % Evaluate the ii box's location in image space 
    %     (if this had not been done before)
    if isempty(Boxes(ii).loc)
        if isempty(Boxes(ii).loc)
            mask=zeros(size(BM));
            mask(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2))=1;
            Boxes(ii).loc=find(mask);
        end
        mask=zeros(size(BM));
        mask(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2))=1;
        Boxes(ii).loc=find(mask);
    end
    
    %loop over the other boxes and find the overlap with ii
    for jj=wh
        %the jj-th box location
        Xb=Boxes(jj).XX;Yb=Boxes(jj).YY;Zb=Boxes(jj).ZZ;
        % Evaluate the ii box's location in image space 
        %     (if this had not been done before)
        if isempty(Boxes(jj).loc)
            mask=zeros(size(BM));
            mask(Xb(1):Xb(2),Yb(1):Yb(2),Zb(1):Zb(2))=1;
            Boxes(jj).loc=find(mask);
        end
        
        % Calculate the overlap between box ii and jj
        [overlap_ii, overlap_jj] = ismember(Boxes(ii).loc,Boxes(jj).loc) ;
        overlap_jj=overlap_jj(overlap_jj>0);
        
        % Calculate the scale factor if they overlap (at least 10% of the voxels)
        if     length( overlap_jj)./length( overlap_ii)>0.1
            Mij=median(Boxes(ii).PD(overlap_ii));
            Mji=median(Boxes(jj).PD(overlap_jj));

            % The scale correction between the two boxes along the overlap region:  boxII= boxJJ*Ratio
            Ratio=median(Boxes(ii).PD(overlap_ii)./ Boxes(jj).PD(overlap_jj));
            
            % Let's check that there is a good overlap, and that after the
            % scale correction the two boxes are similar in the overlap
            % region
            
            % "err" is the median of the sqr error between the two boxes in
            % the overlapregion after the scale correction
            
            %err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii) ./ Boxes(ii).PD(overlap_ii)  ).^2 );
            %err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)).^2 ./ Boxes(ii).PD(overlap_ii)   );
err= median(abs(Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)) ./ Boxes(ii).PD(overlap_ii)  );
            %err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)   ).^2 );
            %err=median(abs(Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii) ./ Boxes(ii).PD(overlap_ii)  ));
            if err<errTres && Ratio>0
               
%                 if Ratio<0.1 || Ratio>3
%                                     keyboard;
%                 end
                % keep the scale Factor between the two boxes
                
                scaleFactor(ii,jj)=1./Ratio;
                scaleFactor(jj,ii)=Ratio;
                
                % Build a linear equation system in which: BoxII - scaleFactor * B0 * JJ =0
                
                % We will aggregate all the boxes like that in LinScaleMat,
                % and solve the linear equation system afterwards
                
%                 LinScaleMat(ii,ii)=1+LinScaleMat(ii,ii);
%                 LinScaleMat(ii,jj)=-Ratio;
            
            end
            
            
        end
        donemask(ii,jj)=1;
        donemask(jj,ii)=1;
    end
    %toc
end
