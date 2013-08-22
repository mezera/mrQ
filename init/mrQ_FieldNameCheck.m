function mrQ=mrQ_FieldNameCheck(mrQ)
% this function cheack competability of the mrQ stracturte. some feilled
% name change and this is a fix to the old starctures;

if (~isfield(mrQ,'SPGR_init_done')) && isfield(mrQ,'SPGR_init') 
    mrQ.SPGR_init_done=mrQ.SPGR_init;
end
if (~isfield(mrQ,'SPGR_coilWeight_done')) && isfield(mrQ,'SPGR_coilWeight')
    mrQ.SPGR_coilWeight_done=mrQ.SPGR_coilWeight;
end
if (~isfield(mrQ,'SPGR_T1fit_done')) && isfield(mrQ,'SPGR_T1fit')
    mrQ.SPGR_T1fit_done=mrQ.SPGR_T1fit;
end
if (~isfield(mrQ,'calM0_done')) && isfield(mrQ,'calM0')
    mrQ.calM0_done=mrQ.calM0;
end
if (~isfield(mrQ,'proclus')) && isfield(mrQ,'proclass') 
    mrQ.proclus=mrQ.proclass;
end

save(mrQ.name,'mrQ')

% in the cases were we re-process and the T1 map file was moved from it original
% location to maps we copy it back so the functions will find it.
if (isfield(mrQ,'spgr_initDir'))
    T1file=fullfile(mrQ.spgr_initDir,'maps/T1_map_lsq.nii.gz');
    if (exist(T1file,'file'))
        
        cmd =(['! cp '    T1file ' ' mrQ.spgr_initDir  '/.']) ;
        eval(cmd);
        
    end
end



