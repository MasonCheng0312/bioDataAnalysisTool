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
        row_With_maxLength["name"] = row_With_maxLength["Transcript Accession"].str.extract(r'(\d+\.\d+|\d+)')
        row_With_maxLength["name"] = row_With_maxLength["name"].astype(float)
        row_With_maxLength.sort_values(by="name",inplace=True)
        row_With_maxLength.drop(columns="name",inplace=True)
        return row_With_maxLength.iloc[[0]]
    
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

def fasta_filter(input_path:str, output_path:str, iso_gene_table:pd.DataFrame) -> None:
    def _checkMode() -> str:
        if iso_gene_table["Gene Type"].iloc[[0]].item() == "PROTEIN_CODING":
            return "p"
        elif iso_gene_table["Gene Type"].iloc[[0]].item() != "PROTEIN_CODING":
            return "n"
        
    def _colName_to_read() -> str:
        if _checkMode() == "p":
            colName = "Protein Accession"
        else:
            colName = "Transcript Accession"
        return colName
    
    def _transformTableToDict() -> dict:
        result = {}
        col = _colName_to_read()
        for _, row in iso_gene_table.iterrows():
            result[row["NCBI GeneID"]] = row[col]
        return result
    
    check_table = _transformTableToDict()
    with open(output_path,"w") as output:
        input_content = SeqIO.parse(input_path,"fasta")
        for record in input_content:
            gene_target = "LOC" + record.description.strip().split("GeneID=")[1].split("]")[0]
            try:
                if check_table[gene_target] == record.name:
                    output.write(">" + record.description + "\n")
                    output.write(f'{record.seq}'+ "\n")
                else:
                    continue
            except:
                continue
    output.close()

if __name__ == "__main__":
    faa_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW3 files/data/Dmel/protein.faa"
    fna_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW3 files/data/Dmel/rna.fna"
    gff_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW3 files/data/Dmel/genomic.gff"
    
    # parse input file
    protein_length_table = parse_protein_seq_length(faa_path)
    protein_length_table.to_csv("/home/cosbi2/py_project/tools/output/protein_length_dmel.csv",index=False)
    transcript_length_table = parse_transcript_seq_length(fna_path)
    # print(transcript_length_table[transcript_length_table["NCBI GeneID"]==30972])
    transcript_length_table.to_csv("/home/cosbi2/py_project/tools/output/trans_length_dmel.csv",index=False)
    protein_transcript_pair = parse_protein_transcript_pair(gff_path)
    protein_transcript_pair.to_csv("/home/cosbi2/py_project/tools/output/pair_dmel.csv",index=False)

    # find isoform
    isoform_candidate = merge_proLen_transLen_to_pairTable(protein_length_table,transcript_length_table,protein_transcript_pair)
    iso_table_for_proteinCoding = find_isoform_for_proteinCoding(isoform_candidate)

    iso_table_for_nonCoding = find_isoform_for_nonCoding(transcript_length_table)
    iso_table_for_nonCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_nonCoding_dmel.csv",index=False)

    iso_table = pd.concat([iso_table_for_nonCoding,iso_table_for_proteinCoding])
    iso_table.to_csv("/home/cosbi2/py_project/tools/output/iso_result_dmel.csv",index=False)
    iso_table_for_proteinCoding.to_csv("/home/cosbi2/py_project/tools/output/iso_result_protein_dmel.csv",index=False) 

    # parse dataset
    dataset_path = "/home/cosbi2/py_project/112_2_bioClassHW/HW3 files/data/Dmel/ncbi_dataset.tsv"
    dsim_table = create_gene_isoform_table(tsv_path=dataset_path,
                                                     protein_iso=iso_table_for_proteinCoding,
                                                     transcript_iso=iso_table_for_nonCoding)
    dsim_table = dsim_table[["NCBI GeneID","Protein Accession","Transcript Accession","Gene Type"]]
    
    # because the geneID col we done doesn't contain "LOC" in front of its num
    dsim_table["NCBI GeneID"] = dsim_table["NCBI GeneID"].apply(add_LOC_to_geneID)

    dsim_table["Protein Accession"] = dsim_table["Protein Accession"].fillna("-")
    dsim_table["Transcript Accession"] = dsim_table["Transcript Accession"].fillna("-")
    
    # from dsim table renew gene fasta file for blast
    blastp_path = "/home/cosbi2/py_project/tools/output/Dmel_for_blastp.fa"
    blastn_path = "/home/cosbi2/py_project/tools/output/Dmel_for_blastn.fa"
    dsim_protein = dsim_table[dsim_table["Gene Type"] == "PROTEIN_CODING"]
    dsim_trans = dsim_table[dsim_table["Gene Type"] != "PROTEIN_CODING"]
    fasta_filter(faa_path,blastp_path,dsim_protein)
    fasta_filter(fna_path,blastn_path,dsim_trans)


    

