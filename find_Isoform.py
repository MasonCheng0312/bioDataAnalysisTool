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
            return result
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

if __name__ == "__main__":
    faa_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/protein.faa"
    fna_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/rna.fna"
    gff_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW2 files/Dsim/genomic.gff"

    protein_length_table = parse_protein_seq_length(faa_path)
    transcript_length_table = parse_transcript_seq_length(fna_path)
    protein_transcript_pair = parse_protein_transcript_pair(gff_path)
    print(protein_transcript_pair[protein_transcript_pair["Protein Accession"] == "XP_039152125.1"])
    print(protein_length_table[protein_length_table["NCBI GeneID"] == 6724725])
    
    # protein_length_table.to_csv("/home/cosbi2/py_project/tools/output/protein_length_table.csv",index=False)
    # transcript_length_table.to_csv("/home/cosbi2/py_project/tools/output/transcript_length_table.csv",index=False)
    # protein_transcript_pair.to_csv("/home/cosbi2/py_project/tools/output/pair_table.csv",index=False)

    isoform_candidate = merge_proLen_transLen_to_pairTable(protein_length_table,transcript_length_table,protein_transcript_pair)
    iso_table_for_proteinCoding = find_isoform_for_proteinCoding(isoform_candidate)

    iso_table_for_nonCoding = find_isoform_for_nonCoding(transcript_length_table)
    iso_table_for_nonCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_nonCoding.csv",index=False)

    iso_table = pd.concat([iso_table_for_nonCoding,iso_table_for_proteinCoding])
    iso_table.to_csv("/home/cosbi2/py_project/tools/output/iso_result.csv",index=False)
    iso_table_for_proteinCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_protein.csv",index=False) 

