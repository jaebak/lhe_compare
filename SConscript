#!/bin/env python3
import os

Import("exportEnv")
envClone = exportEnv.Clone()

source_directories = next(os.walk('src'))[1]
if len(source_directories) == 0:
  source_directories = ['./']

envClone.Append(CPPPATH = ['#/inc'])

# Make binaries for every directory
for source_directory in source_directories: # source_directories = 'core'
  library_objects = []
  # Build objects
  for lib_file in Glob("src/"+source_directory+"*.cpp"): 
    source_object = envClone.Object(lib_file)
    library_objects.append(source_object)
  # Make program
  for source_file in Glob("src/"+source_directory+"/*.cxx"):
    # Build object
    source_object = envClone.Object(source_file)
    # Make binary
    program = envClone.Program('#/kernel/'+envClone['kernel']+'/run/'+source_directory+'/${SOURCE.filebase}.exe', [source_object]+library_objects)
    # Make script that links to binary
    source_basename = os.path.splitext(os.path.basename(str(source_file)))[0]
    source_script_path = 'run/'+source_directory+'/'+source_basename+'.exe'
    source_script_path_mark = "#/"+source_script_path
    envClone.Command(source_script_path_mark, program, './scripts/make_run_scripts.py '+source_script_path)
