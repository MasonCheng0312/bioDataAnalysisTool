import pandas as pd
from Bio import SeqIO

def parse_protein_seq_length(faa_path:str) -> pd.DataFrame:
    file_content = SeqIO.parse(faa_path,"fasta")
    data_of_protein_length_table = []
    for element in file_content:
        # print(element.description)
        # use this to check what it contain in description.
        record = {"NCBI GeneID":element.description.strip().split("GeneID=")[1].split("]")[0],
                  "Protein Accession":element.name,
                  "Protein_length":len(element.seq),}
        data_of_protein_length_table.append(record)
    result = pd.DataFrame(data_of_protein_length_table)
    result.drop_duplicates(keep="first", inplace=True)
    return result

def parse_transcript_seq_length(fna_path:str) -> pd.DataFrame:
    file_content = SeqIO.parse(fna_path,"fasta")
    data_of_transcript_length_table = []
    for element in file_content:
        # print(element.description)
        # use this to check what it contain in description.
        record = {"NCBI GeneID":element.description.strip().split("GeneID=")[1].split("]")[0],
                  "Transcript Accession":element.name,
                  "Transcript_length":len(element.seq),}
        data_of_transcript_length_table.append(record)

    result = pd.DataFrame(data_of_transcript_length_table)
    result.drop_duplicates(keep="first", inplace=True)
    return result

def parse_protein_transcript_pair(gff_path:str) -> pd.DataFrame:
    with open(gff_path,"r") as file:
        content = file.readlines()
        data_of_protein_transcript_pair = []
        for record in content:
            if record.startswith("#"):
                continue
            elif record.split("\t")[2] == "CDS":
                description = record.split("\t")[8]
                data_row = {"Protein Accession":description.split(";")[0].split("ID=cds-")[1],
                            "Transcript Accession":description.split(";")[1].split("Parent=rna-")[1],}
                data_of_protein_transcript_pair.append(data_row)
    result = pd.DataFrame(data_of_protein_transcript_pair) 
    result.drop_duplicates(keep="first", inplace=True) 
    return result

def merge_proLen_transLen_to_pairTable(proLengthTable:pd.DataFrame,transLengthTable:pd.DataFrame,pairTable:pd.DataFrame) -> pd.DataFrame:
    trans = transLengthTable.drop(columns="NCBI GeneID")
    pair_to_protein = pd.merge(proLengthTable, pairTable,on="Protein Accession",how="left")
    trans_to_pair = pd.merge(pair_to_protein,trans,on="Transcript Accession",how="left")
    result = trans_to_pair.reindex(columns=["NCBI GeneID","Protein_length","Protein Accession","Transcript Accession","Transcript_length"])
    return result

def find_isoform_for_proteinCoding(candidateTable:pd.DataFrame) -> pd.DataFrame:
    def _filter(group):
        max_protein_length = group["Protein_length"].max()
        max_protein_row = group[group['Protein_length'] == max_protein_length]
        if len(max_protein_row) == 1:
            return max_protein_row
        else:
            max_trans_length = max_protein_row["Transcript_length"].max()
            result = max_protein_row[max_protein_row['Transcript_length'] == max_trans_length]
            # print(result.iloc[[0]])  
            return result.iloc[[0]]
    iso_table = candidateTable.groupby("NCBI GeneID").apply(_filter).reset_index(drop=True)
    return iso_table

def find_isoform_for_nonCoding(transLenTable:pd.DataFrame) -> pd.DataFrame:
    def _expandTransLengthTable(transLenTable:pd.DataFrame) -> pd.DataFrame:
        transLenTable["Protein Accession"] = "-"
        transLenTable["Protein_length"] = "-"
        return transLenTable

    def _filter(group):
        max_value_of_length = group["Transcript_length"].max()
        row_With_maxLength = group[group["Transcript_length"] == max_value_of_length]
        return row_With_maxLength
    
    trans_candidate_table = _expandTransLengthTable(transLenTable)
    iso_table = trans_candidate_table.groupby("NCBI GeneID").apply(_filter).reset_index(drop=True)
    return iso_table

def create_gene_isoform_table(tsv_path:str, protein_iso:pd.DataFrame, transcript_iso:pd.DataFrame) -> pd.DataFrame:
    def _parse_dataset_tsvFile() -> tuple[pd.DataFrame,pd.DataFrame]:
    # ncbi gene, coding Type
        dataset = pd.read_csv(tsv_path,sep="\t")
        return dataset
    

    data_set= _parse_dataset_tsvFile()

    protein_part = data_set[data_set["Gene Type"] == "PROTEIN_CODING"]
    transcript_part = data_set[data_set["Gene Type"] != "PROTEIN_CODING"]

    protein_iso = protein_iso[["NCBI GeneID","Protein Accession","Transcript Accession"]]
    transcript_iso = transcript_iso[["NCBI GeneID","Protein Accession","Transcript Accession"]]
    protein_part["NCBI GeneID"] = protein_part["NCBI GeneID"].astype('str')
    transcript_part["NCBI GeneID"] = transcript_part["NCBI GeneID"].astype('str')
    protein_result = pd.merge(protein_part,protein_iso,on="NCBI GeneID",how="left")
    trans_result = pd.merge(transcript_part,transcript_iso,on="NCBI GeneID",how="left")

    gene_isoform = pd.concat([protein_result,trans_result])
    return gene_isoform

def add_LOC_to_geneID(name:str) -> str:
    return "LOC"+name

if __name__ == "__main__":
    faa_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/protein.faa"
    fna_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/rna.fna"
    gff_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/genomic.gff"

    protein_length_table = parse_protein_seq_length(faa_path)
    transcript_length_table = parse_transcript_seq_length(fna_path)
    protein_transcript_pair = parse_protein_transcript_pair(gff_path)

    isoform_candidate = merge_proLen_transLen_to_pairTable(protein_length_table,transcript_length_table,protein_transcript_pair)
    iso_table_for_proteinCoding = find_isoform_for_proteinCoding(isoform_candidate)

    iso_table_for_nonCoding = find_isoform_for_nonCoding(transcript_length_table)
    iso_table_for_nonCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_nonCoding.csv",index=False)

    iso_table = pd.concat([iso_table_for_nonCoding,iso_table_for_proteinCoding])
    iso_table.to_csv("/home/cosbi2/py_project/tools/output/iso_result.csv",index=False)
    iso_table_for_proteinCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_protein.csv",index=False) 

    dataset_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/ncbi_dataset.tsv"
    dataset_with_isoform = create_gene_isoform_table(tsv_path=dataset_path,
                                                     protein_iso=iso_table_for_proteinCoding,
                                                     transcript_iso=iso_table_for_nonCoding)
    df = pd.read_csv("/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/output/Dsim/Dsim_table.csv")
    df = df[["NCBI_GeneID","Protein_Accession","Transcript_Accession"]]

    df.rename(columns={"NCBI_GeneID":"NCBI GeneID",
                       "Protein_Accession":"Protein Accession",
                       "Transcript_Accession":"Transcript Accession"},inplace=True)
    dataset_with_isoform = dataset_with_isoform[["NCBI GeneID","Protein Accession","Transcript Accession"]]
    dataset_with_isoform["NCBI GeneID"] = dataset_with_isoform["NCBI GeneID"].apply(add_LOC_to_geneID)
    dataset_with_isoform["Protein Accession"] = dataset_with_isoform["Protein Accession"].fillna("-")
    dataset_with_isoform["Transcript Accession"] = dataset_with_isoform["Transcript Accession"].fillna("-")
    compare = pd.concat([df, dataset_with_isoform]).drop_duplicates(keep= False)
    compare.sort_values(by="NCBI GeneID",inplace=True)
    print(compare)
