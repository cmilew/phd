import numpy as np
import matplotlib.pyplot as plt
import itk
import os,sys,glob


name_file = "*mhd"
path = "./run."
name=[]
folder_list = glob.glob(path + '*')
print(folder_list)
name_idx = 0
for i in range(len(folder_list)):

    file_list = glob.glob( folder_list[i] +"/output**/" + name_file,recursive=True)
    output = "./"
    file_uncertainty = []
    file_edep = []
    file_dose = []
    file_squared = []
    print(len(file_list))
    for i in range(len(file_list)):
        if "uncertainty" in file_list[i]:
            file_uncertainty.append(file_list[i])
        elif "dose" in file_list[i]:
            file_dose.append(file_list[i])
        elif "Squared" in file_list[i]:
            file_squared.append(file_list[i])
        elif "edep" in file_list[i] :
            file_edep.append(file_list[i])

    if len(file_dose) !=0 :
        for i in range(len(file_dose)):
            print(i, file_dose[i])
            if i ==0 :
                name_dose = file_dose[i].split("/")[-1]
            img = itk.imread(file_dose[i])
            array = itk.GetArrayFromImage(img)
            if i ==0 :
                array_tot = np.zeros(array.shape)
            array_tot += array
        img_tot = itk.GetImageFromArray(array_tot)
        img_tot.CopyInformation(img)
        itk.imwrite(img_tot,output + name_dose)

    if len(file_edep) != 0:
        for i in range(len(file_edep)):
            print(i)
            img = itk.imread(file_edep[i])
            array = itk.GetArrayFromImage(img)
            if i == 0:
                name_edep = file_edep[i].split("/")[-1]
            if i ==0 :
                array_tot = np.zeros(array.shape)
            array_tot += array
        img_tot = itk.GetImageFromArray(array_tot)
        img_tot.CopyInformation(img)
        itk.imwrite(img_tot,output + name_edep)





    if len(file_uncertainty) !=0 :
        for i in range(len(file_uncertainty)):
            print(i, file_uncertainty[i])
            err_r_img = itk.imread(file_uncertainty[i])
            err_r_array = itk.GetArrayFromImage(err_r_img)
            if i ==0 :
                name_uncertainty = file_uncertainty[i].split("/")[-1]
            if i==0 :
                err_r_array_tot = np.zeros(err_r_array.shape)
            err_r_array_tot +=err_r_array**2
        err_r_array_tot = np.sqrt(err_r_array_tot)
        err_r_array_tot = err_r_array_tot*1/len(file_uncertainty)

        img_err_r_tot = itk.GetImageFromArray(err_r_array_tot)
        img_err_r_tot.CopyInformation(err_r_img)
        itk.imwrite(img_err_r_tot,output + name_uncertainty)
    name_idx+=1


