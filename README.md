Compares lhe files

Steps

1. Produce LHE files and copy.
  - Untar gridpack
  - `./runcmsgrid.sh 1000 0 8 #1000 is number of events, 0 is random seed, 8 is number of cpus`
  - `scp source_path/cmsgrid_final.lhe target_path/xxx.lhe`

2. source `set_env.sh.mac`

3. scons

4. `./run/compare_lhe.exe 
    -a lhe_files/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.lhe 
    -b lhe_files/ZGTo2LG_2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.lhe 
    --tag_a run2 --tag_b run3 --o run2_run3_ZGToLLG`

5. `open run2_run3_ZGToLLG.pd`
