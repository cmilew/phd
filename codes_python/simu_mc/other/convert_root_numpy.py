import uproot
import glob


path = "./rootfile."
folder_list = glob.glob(path + "*")
name_file = "/*.root"
print(folder_list)
for i in range(len(folder_list)):
    file_path = glob.glob(folder_list[i] + name_file)
    print(f"Num entries for {name_file} :")
    with uproot.open(file_path) as file:
        print(
            "Entries " + file["hits_central_strip;1"]["TotalEnergyDeposit"].num_entries
        )
        print(
            "Deposit " + len(file["hits_central_strip;1"]["TotalEnergyDeposit"].array())
        )
