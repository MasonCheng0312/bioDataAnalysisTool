import pandas as pd


PATH = "/home/cosbi2/orchid/output/Molan/Molan_table.csv"
# table.csv路徑

# ***********將geneID欄位中的"cymsin_Mol029631"的前綴"cymsin_"去除**********
# file_data = pd.read_csv(PATH)
# print(file_data)
# file_data["Gene ID"] = file_data["Gene ID"].str.split("_").str[1]
# print(file_data)
# file_data.to_csv(PATH, index=False)

# ***********將geneID欄位中的"cymsin_Mol029631"的前綴"cymsin_"加回去**********
# file_data = pd.read_csv(PATH)
# file_data["Gene ID"] = "Cymsin_" + file_data["Gene ID"].astype(str)
# file_data["Gene Location"] = "Cymsin_" + file_data["Gene Location"].astype(str)
# file_data["Gene Expression"] = "Cymsin_" + file_data["Gene Expression"].astype(str)
# file_data.to_csv(PATH, index=False)

# ***********刪除多餘欄位***********
# file_data = pd.read_csv(PATH)
# delete_name_of_cols = ["Swissprot_x","Swissprot_y","TrEMBL_x","TrEMBL_y","TrEMBL","GO Term_x"]
# for col in delete_name_of_cols:
#     file_data = file_data.drop(col, axis=1)
# file_data.to_csv(PATH, index=False)

# ***********將sqName欄位中的"cymsin_Mol029631"的前綴"cymsin_"去除**********
# file_data = pd.read_csv(PATH.replace("table", "pep"))
# file_data["sqName"] = file_data["sqName"].str.split("_").str[1]
# print(file_data)
# file_data.to_csv(PATH.replace("table", "pep"), index=False)

# ***********將libraryName欄位重新命名*******************
file_data = pd.read_csv(PATH)
file_data["Library Name"] = "Csinense"
file_data.to_csv(PATH, index=False)
