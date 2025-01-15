import os,glob,shutil
path = './'

run_list = sorted(glob.glob(path + "run.*"))

for run in run_list :
    folder_name = path + "rootfile." +  run.split(".")[-1]
    if os.path.isdir(folder_name):
        shutil.rmtree(folder_name)
    os.mkdir(folder_name)
    rootfile_list = sorted(glob.glob(path + run + "/output**/*.root"))
    if len(rootfile_list) >0 :
        rootfilename = rootfile_list[0].split("/")[-1]
        for i,rootfile in enumerate(rootfile_list):
            if i == 0 :
                cmd_line = "hadd " + folder_name + "/" + rootfilename +" " + rootfile
            else:
                cmd_line += ' ' + rootfile
        os.system(cmd_line)

