function mrQ=mrQ_FieldNameCheck(mrQ)
% this function cheack competability of the mrQ stracturte. some feilled
% name change and this is a fix to the old starctures;

if (~isfield(mrQ,'SPGR_init_done') && sfield(mrQ,'SPGR_init')) 
    mrQ.SPGR_init_done=mrQ.SPGR_init;
end
if (~isfield(mrQ,'SPGR_coilWeight_done') && sfield(mrQ,'SPGR_coilWeight')) 
    mrQ.SPGR_coilWeight_done=mrQ.SPGR_coilWeight;
end
if (~isfield(mrQ,'SPGR_T1fit_done') && sfield(mrQ,'SPGR_T1fit')) 
    mrQ.SPGR_T1fit_done=mrQ.SPGR_T1fit;
end
if (~isfield(mrQ,'calM0_done') && sfield(mrQ,'calM0')) 
    mrQ.calM0_done=mrQ.calM0;
end



