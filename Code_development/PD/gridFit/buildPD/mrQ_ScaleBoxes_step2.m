function [Boxes, scaleFactor]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt)

%mrQ_ScaleBoxes(Boxes,BoxesToUse)
                kk=0;

%LinScaleMat=zeros(32*length(Boxes), length(Boxes));
%LinScaleMat=zeros(length(Boxes), length(Boxes));
scaleFactor=zeros(length(Boxes), length(Boxes));
donemask=ones(length(Boxes), length(Boxes));
donemask(BoxesToUse,BoxesToUse)=0;
%Ref=randperm(BoxesToUse);

BM=readFileNifti(opt.BMfile);
BM=BM.data;
Boxes(1).loc=[];
for ii=BoxesToUse %loop over boxes
    tic
    ii
    donemask(ii,ii)=1;
    % what are the candidate boxes to be overlap with ii box
    wh=find(donemask(ii,:)==0);
    % the ii box location
    XX=Boxes(ii).XX;YY=Boxes(ii).YY;ZZ=Boxes(ii).ZZ;
    
    %evaluate the box ii location in image space if this was not done before
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
    %loop over the other boxes and find the overlap ii
    for jj=wh
        %the jj box location
        Xb=Boxes(jj).XX;Yb=Boxes(jj).YY;Zb=Boxes(jj).ZZ;
        %avaluate the box ii location in image space if this was not done before
        if isempty(Boxes(jj).loc)
            mask=zeros(size(BM));
            mask(Xb(1):Xb(2),Yb(1):Yb(2),Zb(1):Zb(2))=1;
            Boxes(jj).loc=find(mask);
        end
        
        % calculate the overlap between box ii and jj
        [overlap_ii overlap_jj] = ismember(Boxes(ii).loc,Boxes(jj).loc) ;
        overlap_jj=overlap_jj(overlap_jj>0);
        %calculate the scale factor if they overlap (at list 10% of the voxels)
        if     length( overlap_jj)./length( overlap_ii)>0.1
            Mij=median(Boxes(ii).PD(overlap_ii));
            Mji=median(Boxes(jj).PD(overlap_jj));

            % the scale correction between the to boxes along the overlap region  boxII= boxJJ*Ratio
            Ratio=median(Boxes(ii).PD(overlap_ii)./ Boxes(jj).PD(overlap_jj));
            
            % let check that there is a good overlap. and after the scale correction
            % then two box are similar in the overlap region
            %err is the median of the sqr error between the two box in the overlap
            %region after the scale correction
            %err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii) ./ Boxes(ii).PD(overlap_ii)  ).^2 );
             %           err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)).^2 ./ Boxes(ii).PD(overlap_ii)   );
err= median(abs(Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)) ./ Boxes(ii).PD(overlap_ii)  );
%                    err=median((Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii)   ).^2 );
%err=median(abs(Boxes(jj).PD(overlap_jj)*Ratio -Boxes(ii).PD(overlap_ii) ./ Boxes(ii).PD(overlap_ii)  ));
            if err<0.01 && Ratio>0
               
%                 if Ratio<0.1 || Ratio>3
%                                     keyboard;
%                 end
                % keep the scale Factor between the two box
                
                scaleFactor(ii,jj)=1./Ratio;
                scaleFactor(jj,ii)=Ratio;
                
                % build a liner eqation system in wich BoxII -scaleFactor*B0xJJ =0
                % we wll agreagate all the boxes like that in LinScaleMat and solve the linear  eqation system
                % after
                
                %LinScaleMat(ii,ii)=Mij+LinScaleMat(ii,ii);
                %LinScaleMat(ii,jj)=-Mji;
            
            end
            
            
        end
        donemask(ii,jj)=1;
        donemask(jj,ii)=1;
    end
    toc
end
