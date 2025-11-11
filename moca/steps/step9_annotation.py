import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import re

def _standardize_chr_name(series):
    """
    Removes 'chr' prefixes from chromosome names and converts to a lowercase string.
    This ensures consistent chromosome formats (e.g., '1' and 'chr1' are treated as '1').
    """
    return series.astype(str).str.lower().str.replace('^chr', '', regex=True)

def _load_gff(filepath):
    """
    Loads a GFF/GTF file, parsing it into a pandas DataFrame.
    It specifically extracts the gene_id from the attributes column for annotation purposes.
    Handles both GFF (key=value) and GTF (key "value") formats.
    """
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
        # --- Robust Gene ID Extraction Logic ---
        gene_id_col = None

        # 1. Try GFF3 format: ID=gene:xxxxx; or Parent=gene:xxxxx; (often used for 'gene' entries)
        # Note: Your original code's GFF-like extraction seems to target 'gene_id="..."' which is more GTF/GFF-like
        # I will keep the original GFF-like/GTF-like attempts but make them more specific.

        # Try GTF format: gene_id "value";
        gene_id_col = df['attributes'].str.extract(r'gene_id\s*\"([^\"]+)\"', expand=False)
        
        # Try GFF-like format: gene_id=value; or ID=value;
        if gene_id_col.isnull().all():
            print("Could not find 'gene_id \"value\"' attribute. Trying GFF3-like 'gene_id=value;'...")
            gene_id_col = df['attributes'].str.extract(r'gene_id=([^;]+)', expand=False)
            
        # Try GFF3 ID=... format (often for 'gene' or 'mRNA' types)
        if gene_id_col.isnull().all():
            print("Could not find 'gene_id' attribute. Trying GFF3 'ID=...' format...")
            # Capture content until the next semi-colon
            gene_id_col = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)

        # Try GFF3 Parent=... format (for features with a gene parent)
        if gene_id_col.isnull().all():
            print("Could not find 'ID=...' attribute. Trying GFF3 'Parent=gene:...' format...")
            # This is less ideal for the 'gene' row itself, but good for children features.
            gene_id_col = df['attributes'].str.extract(r'Parent=gene:([^;]+)', expand=False)
            
        # --- Final Assignment ---
        df['gene_id'] = gene_id_col
        df['gene_id'] = df['gene_id'].fillna('N/A')
        
        if (df['gene_id'] == 'N/A').all():
             print("\n🚨 **Critical Warning:** Gene ID extraction failed for ALL rows. The pipeline will likely fail later.")
             
        return df
    except Exception as e:
        print(f"Error loading or parsing GFF/GTF file {filepath}: {e}")
        return None

def save_filtered_results(df, output_dir, filter_name):
    """Saves a filtered DataFrame to specifically named CSV and BED files."""
    if df.empty:
        print(f"Warning: No occurrences remained for filter '{filter_name}'. No output files will be generated.")
        return

    # Handle potential column name conflict from merge.
    if 'score_x' in df.columns:
        df = df.rename(columns={'score_x': 'score'}, inplace=False)

    print(f"\nTotal occurrences after '{filter_name}' filter: {len(df)}")
    print(f"Generating output files for '{filter_name}'...")

    # Ensure all required columns are present before saving
    final_cols = ['chr', 'gene_id', 'gene_start', 'gene_end', 'gene_strand', 'motif', 'gen_mstart', 'gen_mend', 'strand', 'score', 'region', 'dist_transc_border']
    final_df_csv = df[[col for col in final_cols if col in df.columns]]

    csv_path = os.path.join(output_dir, f"annotated_motifs_{filter_name}.csv")
    final_df_csv.to_csv(csv_path, index=False, float_format='%.3f')
    print(f"Annotated CSV saved to {csv_path}")

    bed_df = df[['chr', 'gen_mstart', 'gen_mend', 'motif', 'score', 'strand']].copy()
    bed_df.rename(columns={'chr': '#chrom', 'gen_mstart': 'chromStart', 'gen_mend': 'chromEnd', 'motif': 'name'}, inplace=True)
    bed_path = os.path.join(output_dir, f"annotated_motifs_{filter_name}.bed")
    bed_df.to_csv(bed_path, sep='\t', index=False, header=True)
    print(f"Annotated BED file saved to {bed_path}")

def run(config, common_settings):
    """
    Main function for the annotation step.
    Filters motif occurrences based on genomic context and positional preferences.
    """
    # --- 1. Configuration & Path Setup ---
    output_dir = config['output_dir']
    projection_dir = common_settings.get('projection', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    ref_gff_file = config.get('reference_gff')
    chunk_size = config.get('chunk_size', 500000)

    # This flank size MUST match the one used in the extraction step
    EXTRACTION_FLANK_SIZE = 1000 
    
    occurrence_files = sorted(glob.glob(os.path.join(projection_dir, 'occurrences_part_*')))
    species_tag = common_settings.get('species_tag', 'unk')
    model_tag = common_settings.get('model_tag', 'm0')
    tss_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TSS_motif_ranges.csv")
    tts_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TTS_motif_ranges.csv")

    # --- Prepare GFF data for merging ---
    print("Loading and preparing GFF data for merging...")
    gene_annot_df = _load_gff(ref_gff_file)
    if gene_annot_df is None: return

    gene_annot_df = gene_annot_df[gene_annot_df['type'] == 'gene'].copy()
    gene_annot_df.rename(columns={'seqid': 'chr', 'start': 'gene_start', 'end': 'gene_end', 'strand': 'gene_strand'}, inplace=True)
    
    # Create the merge key based on genomic coordinates, matching the occurrence file 'loc' column
    chr_std = _standardize_chr_name(gene_annot_df['chr'])
    start_flank = (gene_annot_df['gene_start'] - EXTRACTION_FLANK_SIZE).astype(str)
    end_flank = (gene_annot_df['gene_end'] + EXTRACTION_FLANK_SIZE).astype(str)
    gene_annot_df['merge_key'] = chr_std + ':' + start_flank + '-' + end_flank
    
    # --- STAGE 1: Process occurrences in chunks to conserve memory ---
    print("\nProcessing occurrence files in chunks...")
    processed_chunks = []
    occ_col_names = ['loc', 'source', 'motif', 'mstart', 'mend', 'score', 'strand', 'pval', 'seq']

    for i, occ_file_part in enumerate(occurrence_files):
        print(f"--> Reading file part {i+1}/{len(occurrence_files)}: {os.path.basename(occ_file_part)}")
        try:
            for occ_chunk in tqdm(pd.read_csv(
                occ_file_part, sep='\t', header=None, low_memory=False,
                names=occ_col_names, chunksize=chunk_size
            ), desc="Processing Chunks"):
                if occ_chunk.empty: continue

                # FIX: Standardize the 'loc' column to match the GFF's merge_key format
                # 1. Split 'loc' into its 'chr' and 'coordinate' parts
                loc_parts = occ_chunk['loc'].str.split(':', n=1, expand=True)
                
                # 2. Standardize the 'chr' part (loc_parts[0]) using the same function as the GFF
                standardized_chr = _standardize_chr_name(loc_parts[0])
                
                # 3. Re-assemble the merge_key
                #    Handle potential missing coordinate part (e.g., if a loc was bad)
                if loc_parts.shape[1] > 1:
                    occ_chunk['merge_key'] = standardized_chr + ':' + loc_parts[1]
                else:
                    # If there's no ':', the loc is invalid; set a key that won't match
                    occ_chunk['merge_key'] = 'invalid_loc_format' 
                # --- End Fix ---
                
                merged_chunk = pd.merge(occ_chunk, gene_annot_df, on='merge_key', how='inner')
                if merged_chunk.empty: continue

                processed_chunks.append(merged_chunk)

        except Exception as e:
            print(f"Error processing file {occ_file_part}: {e}")

    if not processed_chunks:
        print("Warning: No matching motifs found after annotation. Check if 'loc' field in occurrence files matches GFF coordinates.")
        return
        
    # --- Assemble final DataFrame from processed chunks ---
    merged_df = pd.concat(processed_chunks, ignore_index=True)
    print(f"\nSuccessfully processed all chunks. Total occurrences to be filtered: {len(merged_df)}")

    # --- STAGE 2: Calculate Genomic Positions and Define Regions ---
    print("\nStage 2: Calculating genomic positions and applying final filters...")
    flank_start_coord = merged_df['loc'].str.split(':').str[1].str.split('-').str[0].astype(int)
    merged_df['gen_mstart'] = flank_start_coord + merged_df['mstart'] - 1
    merged_df['gen_mend'] = flank_start_coord + merged_df['mend'] - 1
    
    # Define regions based on absolute genomic coordinates
    is_upstream = merged_df['gen_mend'] < merged_df['gene_start']
    is_downstream = merged_df['gen_mstart'] > merged_df['gene_end']
    
    conditions = [is_upstream, is_downstream]
    choices = ['upstream', 'downstream']
    merged_df['region'] = np.select(conditions, choices, default='intragenic')

    # Calculate distance to the nearest transcription border (TSS or TTS)
    dist_choices = [
        merged_df['gene_start'] - merged_df['gen_mend'], # Upstream distance
        merged_df['gen_mstart'] - merged_df['gene_end']   # Downstream distance
    ]
    merged_df['dist_transc_border'] = np.select(conditions, dist_choices, default=0)

    # Correct regions for genes on the minus strand
    is_minus_strand = merged_df['gene_strand'] == '-'
    upstream_on_minus = (merged_df['region'] == 'upstream') & is_minus_strand
    downstream_on_minus = (merged_df['region'] == 'downstream') & is_minus_strand
    merged_df.loc[upstream_on_minus, 'region'] = 'downstream'
    merged_df.loc[downstream_on_minus, 'region'] = 'upstream'
    
    # --- STAGE 3: Apply Positional Preference Filters ---
    try:
        tss_motifs = pd.read_csv(tss_ranges_file)
        tts_motifs = pd.read_csv(tts_ranges_file)
    except FileNotFoundError:
        print(f"Error: Could not read ranging files. Aborting as positional filtering cannot be performed.")
        return

    merged_df['epm'] = merged_df['motif'].str.extract(r'(epm_.+?_p\d+m\d+)', expand=False)
    
    # Separate data by region
    upstream_df = merged_df[merged_df['region'] == 'upstream'].copy()
    downstream_df = merged_df[merged_df['region'] == 'downstream'].copy()
    intragenic_df = merged_df[merged_df['region'] == 'intragenic'].copy() # KEEP intragenic
    
    print(f"Found {len(upstream_df)} upstream, {len(downstream_df)} downstream, and {len(intragenic_df)} intragenic motifs.")

    # Merge with positional preference data
    up_merged = pd.merge(upstream_df, tss_motifs, on='epm', how='left')
    down_merged = pd.merge(downstream_df, tts_motifs, on='epm', how='left')
    
    up_merged.dropna(subset=['min', 'max', 'q10', 'q90'], inplace=True)
    down_merged.dropna(subset=['min', 'max', 'q10', 'q90'], inplace=True)

    # --- Apply 'minmax' filter and save ---
    print("\n--- Applying 'minmax' filter ---")
    up_filtered_minmax = up_merged[up_merged['dist_transc_border'].between(up_merged['min'], up_merged['max'])]
    down_filtered_minmax = down_merged[down_merged['dist_transc_border'].between(down_merged['min'], down_merged['max'])]
    # Combine filtered upstream/downstream motifs WITH ALL intragenic motifs
    final_df_minmax = pd.concat([up_filtered_minmax, down_filtered_minmax, intragenic_df], ignore_index=True)
    save_filtered_results(final_df_minmax, output_dir, "minmax")

    # --- Apply 'q1q9' filter and save ---
    print("\n--- Applying 'q1q9' filter ---")
    up_filtered_q1q9 = up_merged[up_merged['dist_transc_border'].between(up_merged['q10'], up_merged['q90'])]
    down_filtered_q1q9 = down_merged[down_merged['dist_transc_border'].between(down_merged['q10'], down_merged['q90'])]
    # Combine filtered upstream/downstream motifs WITH ALL intragenic motifs
    final_df_q1q9 = pd.concat([up_filtered_q1q9, down_filtered_q1q9, intragenic_df], ignore_index=True)
    save_filtered_results(final_df_q1q9, output_dir, "q1q9")

    print("\nAnnotation step successfully completed.")
