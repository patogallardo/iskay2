all: params.json
	iskay2_mk_paramfiles_sweep_lum.py
	iskay2_mk_sbatch_files.py
	sbatch sbatch_DR6_150GHz_C21cat_lum_gt_4p3e10.sh
clean:
	$(RM) results/*.hdf
