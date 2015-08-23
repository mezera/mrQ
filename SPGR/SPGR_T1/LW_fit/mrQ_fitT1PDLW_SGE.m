function mrQ_fitT1PDLW_SGE(IDnum,jumpindex,jobindex)
%
% function mrQ_fitT1PDLW_SGE(opt,jumpindex,jobindex)
%
% Perform the T1 and PD fitting, using the SunGrid
%
%   ~INPUTS~
%        IDnum: A unique numerical identifier for each subject, it is
%                    created by taking the first eight digits of the
%                    temporary file name
%    jumpIndex: The number of voxels to send for each job per machine on
%                    the grid. If SunGrid is not being used (that is, only
%                    one machine is available), then the jumpIndex will
%                    just be all the voxels.
%     jobIndex: The number of the job.
%
% Saves: 'res','resnorm','st','ed'
%    TO: [opt.outDir opt.name '_' num2str(st) '_' num2str(ed)]
%
% See Also:
%   mrQ_fitT1M0.m, mrQ_fitT1PD_LSQ.m
%
% Authors: Aviv Mezer and Nikola Stikov date 24.02.2014



%% I. Get the opt structure using the mrQ ID
mrQpath= mrQ_getPath(IDnum);
load(mrQpath);
load(mrQ.LWoptname);

%% II. Set the maximum number of computational threads available to Matlab
%maxnumcompthreads(1)

j  = 0;
st = 1 +(jobindex-1)*jumpindex;
ed = st+jumpindex-1;

if ed>length(opt.wh)
    ed = length(opt.wh);
end

% parfor can be used in here
for i= st:ed,
    j=j+1;
    
    if (~isempty(find(isnan(opt.s(i,:)))) | ~isempty(find(opt.s(i,:)==0)) |   ~isempty(find(isinf(opt.s(i,:)))))
        res(1:4,j)=nan;
        
    else
        TR=opt.tr;
        fa=opt.flipAngles/180*pi*opt.B1(i);
        y = (opt.s(i,:))./sin(fa);
        x = (opt.s(i,:))./tan(fa);
        Vals = polyfit(x, y, 1);
        slopeBiased = Vals(1);
        result = abs(-TR./log(slopeBiased));
        
        %slopeBiased=pi^(-TR/linearT1)
        %
        t1Biased= result;
        
        %% PD
        %biased
        pdBias=opt.s(i,:)./((1-exp(-TR./t1Biased)).*sin(fa)./(1-exp(-TR./t1Biased).*cos(fa)));
        pdBias=mean(pdBias)./opt.Gain(i);
        
        weights = (sin(fa)./(1 - slopeBiased.*cos(fa))).^2;
        weights(isinf(weights)) = 0; %remove points with infinite weight
        
        Vals2 = polyfitweighted(x, y, 1, weights);
        slopeUnbiased = Vals2(1);
        t1Unbiased = abs(-TR./log(slopeUnbiased));
        
        %unbiased
        
        pdUnBias=opt.s(i,:)./((1-exp(-TR./t1Unbiased)).*sin(fa)./(1-exp(-TR./t1Unbiased).*cos(fa)));
%         pdUnBias1=mean(pdUnBias)./opt.Gain(i);
        pdUnBias=(sum(weights.*pdUnBias))./(opt.Gain(i).*sum(weights));
        res(1,j)=t1Unbiased;
        res(2,j)=pdUnBias;
        res(3,j)=t1Biased;
        res(4,j)=pdBias;
        
        
        
    end
end
List={'t1Unbiased' 'pdUnBias','t1biased' 'pdBias'};
%%

name = [opt.outDir opt.name '_' num2str(st) '_' num2str(ed)];
save(name,'res','st','ed','List')
