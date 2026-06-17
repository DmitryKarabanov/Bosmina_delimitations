import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# 1. Region setup (order matters: it defines columns in the TRAITS block)
country_to_region = {
    'Canada': 'North_America', 'USA': 'North_America',
    'Mexico': 'Central_America', 'Guatemala': 'Central_America', 'Panama': 'Central_America',
    'Colombia': 'South_America', 'Ecuador': 'South_America', 'Peru': 'South_America',
    'Brazil': 'South_America', 'Argentina': 'South_America', 'Chile': 'South_America',
    'Spain': 'Europe', 'Czechia': 'Europe', 'Finland': 'Europe', 'Italy': 'Europe',
    'Israel': 'Asia_Middle_East', 'Japan': 'Asia_Middle_East', 'China': 'Asia_Middle_East', 'Russia': 'Asia_Middle_East',
    'Australia': 'Australia'
}
trait_labels = ['North_America', 'Central_America', 'South_America', 'Europe', 'Asia_Middle_East', 'Australia']
ntraits = len(trait_labels)

# 2. Load CSV
print("Reading CSV...")
df = pd.read_csv('_resources.csv', sep=';')

# In your CSV, the ID is in the 'seq' column
id_col = 'seq' if 'seq' in df.columns else df.columns[1]
df[id_col] = df[id_col].astype(str).str.strip()
df['County'] = df['County'].astype(str).str.strip()

seq_to_region = {}
for _, row in df.iterrows():
    seq_id = str(row[id_col]).strip()
    country = str(row['County']).strip()
    
    # Skip technical rows and empty values
    if seq_id in ['nan', ''] or 'GenBank' in seq_id or 'acc' in seq_id.lower():
        continue
        
    if country in country_to_region:
        seq_to_region[seq_id] = country_to_region[country]

print(f"Successfully mapped to regions: {len(seq_to_region)} sequences from CSV.")

# 3. Read FASTA and collapse into haplotypes
print("Processing FASTA and collapsing identical sequences...")
haplotypes = defaultdict(list) # {sequence_string: [list_of_seq_ids]}
unmatched_fasta = []

for record in SeqIO.parse('_bosmina_no_root.fasta.txt', 'fasta'):
    # Take the FULL ID, as it completely matches the 'seq' column in CSV!
    clean_id = record.id.strip()
    seq_str = str(record.seq).upper().replace('\n', '').replace('\r', '').replace(' ', '').replace('-', '')
    
    if clean_id in seq_to_region:
        haplotypes[seq_str].append(clean_id)
    else:
        unmatched_fasta.append(clean_id)

if unmatched_fasta:
    print(f"Warning: {len(unmatched_fasta)} sequences from FASTA were not found in CSV (or have no region).")
    print(f"   Examples: {unmatched_fasta[:3]}")

if not haplotypes:
    print("Critical error: No matches found between FASTA and CSV!")
else:
    print(f"Found matches: {sum(len(v) for v in haplotypes.values())} sequences collapsed into {len(haplotypes)} unique haplotypes.")

# Fix the order of haplotypes for strict correspondence between DATA and TRAITS blocks
hap_data = []
for i, (seq_str, ids) in enumerate(haplotypes.items(), 1):
    hap_data.append((f"H{i:03d}", seq_str, ids))

ntax = len(hap_data)

if ntax == 0:
    print("Haplotype list is empty. NEXUS file not created.")
else:
    seq_lengths = [len(s) for _, s, _ in hap_data]
    nchar = max(seq_lengths)
    
    if len(set(seq_lengths)) > 1:
        print(f"Warning: Sequence length varies ({min(seq_lengths)}-{max(seq_lengths)} bp).")
        print("   Shorter ones will be padded with '?' for correct PopART processing.")

    # 4. Generate NEXUS
    print("Generating popart_input.nex...")
    nexus_lines = []
    nexus_lines.append("#NEXUS")
    nexus_lines.append("")
    nexus_lines.append("BEGIN DATA;")
    nexus_lines.append(f"  Dimensions NTAX={ntax} NCHAR={nchar};")
    nexus_lines.append("  Format DataType=DNA Gap=- Missing=?;")
    nexus_lines.append("  Matrix")
    
    trait_matrix_rows = []
    
    for hap_name, seq_str, ids in hap_data:
        padded_seq = seq_str.ljust(nchar, '?')
        nexus_lines.append(f"  {hap_name}  {padded_seq}")
        
        region_counts = {r: 0 for r in trait_labels}
        for sid in ids:
            reg = seq_to_region.get(sid)
            if reg in region_counts:
                region_counts[reg] += 1
        
        trait_row = ", ".join(str(region_counts[r]) for r in trait_labels)
        trait_matrix_rows.append(f"  {hap_name} {trait_row}")
    
    nexus_lines.append("  ;")
    nexus_lines.append("END;")
    nexus_lines.append("")
    nexus_lines.append("BEGIN TRAITS;")
    nexus_lines.append(f"  Dimensions NTRAITS={ntraits};")
    nexus_lines.append("  Format labels=yes missing=? separator=Comma;")
    nexus_lines.append(f"  TraitLabels {' '.join(trait_labels)};")
    nexus_lines.append("  Matrix")
    nexus_lines.extend(trait_matrix_rows)
    nexus_lines.append("  ;")
    nexus_lines.append("END;")
    
    with open("popart_input.nex", "w", encoding="utf-8") as f:
        f.write("\n".join(nexus_lines) + "\n")
    
    print(f"Done! Created 'popart_input.nex' with {ntax} unique haplotypes.")