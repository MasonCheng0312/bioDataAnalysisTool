import pandas as pd
import numpy as np
import re

class blastParser():
    def __init__(self, path1:str, path2:str, mainSpecies:str, minorSpecies:str) -> None:
        # notes the name of main species and minor species.
        self.main = mainSpecies
        self.minor = minorSpecies

        # read csv file to dataframe and give it column name.
        self.colName1 = [mainSpecies, "E value", "Score", "Query coverage", "identity", minorSpecies, "UK1", "UK2"]
        self.file1 = pd.read_csv(path1,sep="\t",names=self.colName1)
        self.colName2 = [minorSpecies, "E value", "Score", "Query coverage", "identity", mainSpecies, "UK1", "UK2"]
        self.file2 = pd.read_csv(path2,sep="\t",names=self.colName2)

        # delete "UK1"&"UK2" column.
        self.colName1 = self.colName1[:6]
        self.colName2 = self.colName2[:6]

    def bestMatchFilter(self, e_value) -> tuple[pd.DataFrame,pd.DataFrame]:
        def _bestHitChooser(data,minor) -> pd.DataFrame:  # for each group data, choose the best hit one(or more).
            row_of_minE_Value = data[data["E value"] == data["E value"].min()]
            if len(row_of_minE_Value) == 1:
                return row_of_minE_Value.values.tolist()
            else:
                row_of_maxScore = row_of_minE_Value[row_of_minE_Value["Score"] == row_of_minE_Value["Score"].max()]
                if len(row_of_maxScore) == 1:
                    return row_of_maxScore.values.tolist()
                else:
                    row_of_maxQuery = row_of_maxScore[row_of_maxScore["Query coverage"] == row_of_maxScore["Query coverage"].max()]
                    if len(row_of_maxQuery) == 1:
                        return row_of_maxQuery.values.tolist()
                    else:
                        row_of_maxIdentity = row_of_maxQuery[row_of_maxQuery["identity"] == row_of_maxQuery["identity"].max()]
                        if len(row_of_maxIdentity) == 1:
                            return row_of_maxIdentity.values.tolist()
                        else:
                            # this part will merge all row with same condition(e-value,query coverage...etc) to one.
                            # name = []
                            # for row in row_of_maxIdentity.values:
                            #     name.append(row[0].split(" ")[-1])
                            # result = row_of_maxIdentity.iloc[0]
                            # result[minor] = ",".join(name)
                            # return result.values.tolist()
                            return row_of_maxIdentity.values.tolist()

                        
        # first, we drop the last 2 col of csv.
        col_to_delete = ["UK1", "UK2"]
        self.file1.drop(col_to_delete, inplace=True, axis=1)
        self.file2.drop(col_to_delete, inplace=True, axis=1)

        # then, we delete the row which e-value bigger than our input limit.
        self.file1 = self.file1[self.file1["E value"] < e_value]
        self.file2 = self.file2[self.file2["E value"] < e_value]

        # finally, we want to get the best hit of blast result
        group1 = self.file1.groupby(self.main)
        result1 = []
        group2 = self.file2.groupby(self.minor)
        result2 = []

        for _, group_data in group1:
            record = _bestHitChooser(group_data,self.minor)
            for element in record:
                result1.append(element[0:6])

        for group_name, group_data in group2:
            record = _bestHitChooser(group_data,self.main)
            for element in record:
                result2.append(element[0:6])

        return pd.DataFrame(result1,columns=self.colName1), pd.DataFrame(result2,columns=self.colName2)

    def handleBlastData(self, main_to_minor:pd.DataFrame, minor_to_main:pd.DataFrame):
        merge_table = pd.merge(main_to_minor,minor_to_main,on=self.minor,how="left")

        # classify the type of merge result.
        merge_table['type'] = np.where(
            (merge_table['Dsim_x'] == merge_table['Dsim_y']) & (~merge_table['Dsim_y'].isnull()), # 條件1：Dsim_x 等於 Dsim_y 且 Dsim_y 不為 NaN
            1,  # 滿足條件1時的值為 1
            np.where(
                merge_table['Dsim_y'].isnull(),  # 條件2：Dsim_y 為 NaN
                3,  # 滿足條件2時的值為 3
                2   # 不滿足條件1和條件2時的值為 2
            )
        )
        merge_table.rename(columns={"Dsim_x":"Dsim","Dsim_y":"hit back"},inplace=True)
        return merge_table
    
def parseGeneTableForBlast(tablePath:str, speciesName:str) -> tuple[pd.DataFrame,pd.DataFrame]:
    geneTable = pd.read_csv(tablePath,sep=",")
    dropCol = ["Symbol", "Description", "Transcripts", "Nomenclature ID", "Nomenclature ID2", "Chromosomes", "Transcript_Genomic_Accession", "Transcript_Genomic_Start", "Transcript_Genomic_Stop", "Orientation", "Proteins", "Synonyms"]
    geneTable.drop(columns=dropCol, inplace=True, axis=1)

    table_for_blast_n = geneTable.drop(columns="Protein_Accession").rename(columns={"Transcript_Accession":speciesName})
    table_for_blast_n = table_for_blast_n[table_for_blast_n[speciesName] != "-"]

    table_for_blast_p = geneTable.drop(columns="Transcript_Accession").rename(columns={"Protein_Accession":speciesName})
    table_for_blast_p = table_for_blast_p[table_for_blast_p[speciesName] != "-"]

    if speciesName == "Dsim":
        table_for_blast_n.rename(columns={"NCBI_GeneID":"Gene name"},inplace=True)
        table_for_blast_p.rename(columns={"NCBI_GeneID":"Gene name"},inplace=True)

    elif speciesName == "Dmel":
        table_for_blast_n.rename(columns={"NCBI_GeneID":"Dmel name"},inplace=True)
        table_for_blast_n.drop(columns="Gene_Type",inplace=True)

        table_for_blast_p.rename(columns={"NCBI_GeneID":"Dmel name"},inplace=True)
        table_for_blast_p.drop(columns="Gene_Type",inplace=True)

    return table_for_blast_n, table_for_blast_p

def mergeTable(blastResult:pd.DataFrame, main_species_table:pd.DataFrame, minor_species_table:pd.DataFrame) -> pd.DataFrame:
    main = "Dsim"
    minor = "Dmel"

    merge_result = pd.merge(blastResult,main_species_table, left_on= main,right_on=main,how="left")
    # print(merge_result)
    main_species_table.drop(columns="Gene_Type",inplace=True)
    main_species_table.rename(columns={"Gene name":"Hit back name"},inplace=True)
    merge_result2 = pd.merge(merge_result,main_species_table, right_on=main, left_on="hit back",how="left")

    merge_result3 = pd.merge(merge_result2,minor_species_table,on=minor,how="left")
    merge_result3.drop(columns="hit back",inplace=True)
    merge_result3.rename(columns={'Dsim_x': 'Dsim', 'E value_x': 'evalue1', 'Score_x': 'score1', 'Query coverage_x': 'qcover1', 'identity_x': 'identity1', 'E value_y': 'evalue2', 'Score_y': 'score2', 'Query coverage_y': 'qcover2', 'identity_y': 'identity2', 'Dsim_y': 'hit_back', 'type': 'Type'}, inplace=True)

    return merge_result3

def parseGeneTableForOrtholog(tablePath:str) -> pd.DataFrame :
    geneTable = pd.read_csv(tablePath,sep=",")
    dropCol = ["Symbol", "Description", "Transcripts", "Nomenclature ID", "Nomenclature ID2", "Chromosomes",
                "Transcript_Genomic_Accession", "Transcript_Genomic_Start", "Transcript_Genomic_Stop", "Orientation",
                "Proteins", "Synonyms", "Protein_Accession", "Transcript_Accession"]
    geneTable.drop(columns=dropCol, inplace=True, axis=1)
    geneTable.rename(columns={"NCBI_GeneID":"Gene name"}, inplace=True)
    
    Dsim_protein = geneTable[geneTable["Gene_Type"] == "PROTEIN_CODING"]
    Dsim_transcript = geneTable[geneTable["Gene_Type"] != "PROTEIN_CODING"]

    return Dsim_protein, Dsim_transcript

def handleCheckTable(blastResult:pd.DataFrame) -> pd.DataFrame:
    def list_elements_to_str(lst):
        return list(map(str, lst))
    
    blastResult.drop(columns=["evalue1", "evalue2", "score1", "score2", "qcover1", "qcover2", "identity1", "identity2", "Dsim", "Dmel", "Gene_Type", "hit_back","Hit back name"],inplace=True)
    # print(blastResult)
    # blastResult.drop_duplicates(keep="first",inplace=True)
    Type1_table = blastResult[blastResult["Type"] == 1]
    # Type1_table.drop_duplicates(keep="first",inplace=True)
    Type1_result = Type1_table.groupby("Gene name").agg({
                                                        #  "Gene name":"first", 
                                                         "Type":"size", 
                                                         "Dmel name":lambda x: x.tolist()}).reset_index()
    # print(Type1_result)
    Type1_result.rename(columns={"Type":"Ortholog_Num"},inplace=True)
    Type1_result["Dmel name"] = Type1_result["Dmel name"].apply(list_elements_to_str)

    Alltype_result = blastResult.groupby("Gene name").agg({"Type":lambda x: x.tolist()}).reset_index()
    result = pd.merge(Alltype_result,Type1_result,on="Gene name",how="left")
    result["Ortholog_Num"].fillna(0.0, inplace=True)

    return result
    
def merge_GeneTable_and_BlastResult(geneTable:pd.DataFrame,blastResult:pd.DataFrame,namelist:list) -> pd.DataFrame:
    def check_gene_name(row): 
    # 如果 'Type' 列的值是 NaN，检查 'Gene name' 是否存在于名单中
        if pd.Series(pd.isna(row['Type'])).any():
            if row['Gene name'] in namelist:
                return '[4]'  # 如果 'Type' 列是 NaN 并且 'Gene name' 列在名单中，则返回 [4]
            else:
                return '-'
        else:
            return row['Type']  # 否则返回 'Type' 列的原始值
        
    result = pd.merge(geneTable,blastResult,on="Gene name",how="left")
    result["Type"] = result.apply(check_gene_name, axis=1)
    result["Ortholog_Num"].fillna(0.0, inplace=True)
    return result

def findGeneNameWithSeq() -> list:
    with open("/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/Dsim_for_blastn.fa","r") as file:
        content = file.readlines()
    gene_names = []
    # 使用正则表达式从每一行中提取基因名称
    for line in content:
        match = re.search(r'LOC(\d+)', line)
        if match:
            gene_names.append(match.group())
    file.close()
    with open("/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/Dsim_for_blastp.fa","r") as file:
        content = file.readlines()
    for line in content:
        match = re.search(r'LOC(\d+)', line)
        if match:
            gene_names.append(match.group())
    file.close()

    # 去除重复的基因名称
    gene_names = list(set(gene_names))
    return gene_names

if __name__ == "__main__":
    main = "Dsim"
    minor = "Dmel"

    # path for blastn
    PATH1 = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/blastn_Dsim_inDB_Dmel.csv"
    PATH2 = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/blastn_Dmel_inDB_Dsim.csv"

    blastN_processor = blastParser(PATH1,PATH2,main,minor)
    bestHit_Dsim_to_Dmel, bestHit_Dmel_to_Dsim = blastN_processor.bestMatchFilter(1e-5)
    blast_n_result = blastN_processor.handleBlastData(bestHit_Dsim_to_Dmel, bestHit_Dmel_to_Dsim)

    # path for blastp
    PATH3 = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/blastp_Dsim_inDB_Dmel.csv"
    PATH4 = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/blastp_Dmel_inDB_Dsim.csv"

    blastP_processor = blastParser(PATH3,PATH4,main,minor)
    bestHit_Dsim_to_Dmel, bestHit_Dmel_to_Dsim = blastP_processor.bestMatchFilter(1e-5)
    blast_p_result = blastP_processor.handleBlastData(bestHit_Dsim_to_Dmel, bestHit_Dmel_to_Dsim)

    # parse Gene Table to delete column we dont need.
    DsimTablePath = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/Dsim_table.csv"
    DsimTable_blastN, DsimTable_blastP = parseGeneTableForBlast(DsimTablePath,"Dsim")

    DmelTablePath = "/home/cosbi2/py_project/112_2_bioClassHW/HW1 files/input/Dmel_table.csv"
    DmelTable_blastN, DmelTable_blastP = parseGeneTableForBlast(DmelTablePath,"Dmel")
    
    check_transcript_table = mergeTable(blast_n_result,DsimTable_blastN,DmelTable_blastN)
    check_protein_table = mergeTable(blast_p_result,DsimTable_blastP,DmelTable_blastP)
 
    # print(check_protein_table)
    # result = check_protein_table.drop_duplicates(keep="first",inplace=False,subset=None)
    # print(result)
    # df = pd.read_csv("112_2_bioClassHW/HW1 files/output/Dsim_check_protein.csv")
    # print(df)
    # df2 = df.drop_duplicates(keep="first",inplace=False,subset=None)

    Dsim_protein, Dsim_transcript = parseGeneTableForOrtholog(DsimTablePath)
    # print(Dsim_Table)
    proteinTable = handleCheckTable(check_protein_table)
    transcriptTable = handleCheckTable(check_transcript_table)                  
    
    totalNamelist = findGeneNameWithSeq()
    protein_result = merge_GeneTable_and_BlastResult(Dsim_protein, proteinTable,totalNamelist)
    transcript_result = merge_GeneTable_and_BlastResult(Dsim_transcript, transcriptTable,totalNamelist)

    Finalresult = pd.concat([protein_result,transcript_result])
    Finalresult.to_csv("/home/cosbi2/py_project/tools/output/result.csv",index=False)
    # protein_result.to_csv("/home/cosbi2/py_project/tools/output/result.csv",index=False)


