function mrQ=mrQ_Path_checks(mrQ)
% mrQ=mrQ_Path_checks(mrQ)
% This function is used to make v.2  and v.2.1 compatible. There are some
% different in where path are saved in mrQ between the version.
% this is a patch that allow the part that after the B1 fit and
% registration to work without a change in the code.
%
%
% (C) A.M. Mezer lab, the Hebrew University of Jerusalem, Israel
% 2016
%%
if ~isfield(mrQ,'spgr_initDir') % this is needed for continuation between mrQ versions
    mrQ.spgr_initDir=mrQ.InitSPGR.spgr_initDir;
    mrQ.xform=mrQ.InitSPGR.xform;
    mrQ.HeadMask=mrQ.LinFit.HeadMask;
    mrQ.BrainMask=mrQ.LinFit.BrainMask;
    mrQ.T1_LFit=mrQ.LinFit.T1_LFit;
    mrQ.M0_LFit=mrQ.LinFit.M0_LFit;
    mrQ.T1_LFit_HM=mrQ.LinFit.T1_LFit_HM;
    mrQ.M0_LFit_HM=mrQ.LinFit.M0_LFit_HM;
    
    
    
end
