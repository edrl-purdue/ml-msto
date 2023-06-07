% SETUP PARALLEL POOL
function setup_parpool()
    pc = parcluster('local'); %process based
	dpc = pc.JobStorageLocation;
	if ~isempty(getenv('SLURM_CPUS_ON_NODE'))
		rmdir(dpc,'s')
		delete(gcp('nocreate'))
		parpool(pc,str2double(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout',480);
	else
		if ~ispc
			rmdir(dpc,'s')
			delete(gcp('nocreate'))
			parpool(pc,feature('numcores'),'IdleTimeout',480);
		end
	end
end