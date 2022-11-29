;Section 8, using IDL
ctm_get_data, datainfo, filename='/home/grifbake/chen3274/bbyrne/runs/v8-02-01/m9/OptData/gctm.gdt.01','IJ-GDE-$', tracer=91003
ctm_read_data, data, datainfo
write_csv,'/home/grifbake/chen3274/bbyrne/runs/v8-02-01/m9/OptData/gdt.csv',data
exit


