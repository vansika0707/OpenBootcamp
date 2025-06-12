import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the TSV file
df = pd.read_csv("eukaryotes.tsv", sep="\t")

# Show basic descriptive statistics for Size (Mb)
print("Descriptive statistics for 'Size (Mb)':")
print(df["Size (Mb)"].describe())

# Filter rows where Species is Homo sapiens
humans = df[df["Species"] == "Homo sapiens"]
print("\nFiltered DataFrame for Homo sapiens:")
print(humans)

# Convert Size (Mb) and Number of genes to numeric, handle non-numeric entries like '-'
df["Size (Mb)"] = pd.to_numeric(df["Size (Mb)"], errors="coerce")
df["Number of genes"] = pd.to_numeric(df["Number of genes"], errors="coerce")

# Compute correlation (will skip rows with NaN)
correlation = df[["Size (Mb)", "Number of genes"]].corr()
print("\nCorrelation between 'Size (Mb)' and 'Number of genes':")
print(correlation)

# Create a DataFrame of small genomes (Size < 5000)
small_genomes = df[df["Size (Mb)"] < 5000]

# Drop rows with NaNs before plotting
small_genomes = small_genomes.dropna(subset=["Size (Mb)", "Number of genes"])

# Plot: Number of genes vs Size (Mb)
plt.figure(figsize=(10, 6))
sns.scatterplot(data=small_genomes, x="Size (Mb)", y="Number of genes")
plt.title("Gene Count vs Genome Size (<5000 Mb)")
plt.xlabel("Genome Size (Mb)")
plt.ylabel("Number of Genes")
plt.tight_layout()
plt.savefig("genome_plot.png")  # Save the plot
