import pandas as pd

def calculate_pheno_stats(subset, gene, total):
    pheno_stats = subset.iloc[:, [3]].value_counts()
    pheno_df = pheno_stats.to_frame().reset_index()
    pheno_df.columns = ["Phenotype", "count"]
    pheno_df["Gene"] = gene
    
    # Handle "Phenotype could not be inferred"
    inferred_row = pheno_df[pheno_df["Phenotype"] == "Phenotype could not be inferred"]
    
    if not inferred_row.empty:
        result = inferred_row.iloc[0]["count"]
        updated_total = total - result
        pheno_df["Phenotype could not be inferred"] = result
        pheno_df.drop(inferred_row.index, inplace=True)
    else:
        pheno_df["Phenotype could not be inferred"] = 0

    pheno_df["Total_count"] = pheno_df["count"]
    pheno_df["Percent"] = (pheno_df["count"] / updated_total) * 100 if 'updated_total' in locals() else (pheno_df["count"] / total) * 100
    
    return pheno_df

def calculate_haplo_stats(subset, gene, total):
    # Combine haplotype columns
    combined_column = pd.concat([subset.iloc[:, 0], subset.iloc[:, 1]], ignore_index=True)
    
    # Create a DataFrame from the combined column
    haplo_df = pd.DataFrame({'combined': combined_column})
    
    # Use groupby to count haplotype occurrences
    stats_df = haplo_df.groupby('combined').size().reset_index(name='count')
    
    stats_df["Gene"] = gene
    stats_df["AF"] = stats_df["count"] / (total * 2)
    stats_df["AF_Percent"] = stats_df["AF"] * 100
    
    return stats_df

# Load data
df = pd.read_csv("V3_combine.csv")
genes = ["ABCG2", "CFTR", "CFTR1", "CFTR2", "CYP2B6", "CYP2C19", 
         "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5", "CYP4F2", 
         "DPYD", "G6PD", "IFNL3", "MTRNR1", "NAT2", 
         "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", 
         "VKORC1"]

total_samples = 4134
freq_output = []
stats_output = []

for gene in genes:
    subset = df.loc[:, df.columns.str.startswith(gene)]
    
    freq_output.append(calculate_pheno_stats(subset, gene, total_samples))
    stats_output.append(calculate_haplo_stats(subset, gene, total_samples))

# Concatenate results and save to CSV files
freq_df = pd.concat(freq_output)
df_stats = pd.concat(stats_output)

freq_df.to_csv("phenotype_frequency.csv", index=False)
df_stats.to_csv("allele_frequency.csv", index=False)
