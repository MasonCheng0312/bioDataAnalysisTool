import pandas as pd
from Bio import SeqIO

input_path = "/home/cosbi2/py_project/input/JL_info.csv"
jlan_data = pd.read_csv(input_path,sep="\t",usecols=['qseqID', 'pident', 'E value', 'bitscore'])
# jlan_data = jlan_data['qseqID', 'pident', 'E value', 'bitscore']
jlan_data.rename(columns={"pident":"Identities","bitscore":"Score","E value":"Expect value","qseqID":"Gene ID"},inplace=True)
print(jlan_data)

origin_table_path = "/home/cosbi2/py_project/input/Censifolium_table.csv"
jlan_col = ['Library Name', 'Gene ID', 'Gene Location', 'Gene Expression',
       'GenBank', 'Accession number(Best hits in the GenBank)', 'Annotation',
       'Species', 'Score', 'Expect value', 'Identities', 'Frame',
       'KEGG Pathway', 'GO Term', 'Interpro', 'Swissprot', 'TrEMBL', 'mature miRNA']
jlan_table = pd.read_csv(origin_table_path, sep=",", index_col=None,header=None)
jlan_table.columns = jlan_col
jlan_table.drop(columns=['Score', 'Expect value', 'Identities'],inplace=True)
# print(jlan_table)

result = pd.merge(jlan_table,jlan_data,on="Gene ID",how="left")
# print(result)
result = result.reindex(columns= jlan_col)
result.to_csv("/home/cosbi2/py_project/tools/output/JL_result.csv",index=False)
data_to_handle = result[result["Expect value"].isna()]
name_list = data_to_handle["Gene ID"].tolist()
# print(name_list)
input_pep = "/home/cosbi2/py_project/input/Censifolium_pep.fa"
outputPATH = "/home/cosbi2/py_project/tools/output/mismatch_pep.fa"
with open(outputPATH,"w") as output:
	input_content = SeqIO.parse(input_pep,"fasta")
	for record in input_content:
		if record.name in name_list:
			output.write(">" + record.name + "\n")
			output.write(f'{record.seq}'+ "\n")
		else:continue
