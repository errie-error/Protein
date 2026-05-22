import argparse

from protein_project.config import ensure_project_dirs, get_processed_path, get_raw_path, load_config
from protein_project.data import clean_clinvar_table, download_file, load_variant_summary, read_fasta_sequence
from protein_project.structure import assign_spatial_clusters, build_saprot_sequences, extract_residue_table, save_saprot_sequences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/tp53.yaml")
    parser.add_argument("--foldseek", default=None)
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)

    clinvar_path = download_file(
        config["urls"]["clinvar_variant_summary"],
        get_raw_path(config, "clinvar_variant_summary"),
        force=args.force_download,
    )
    fasta_path = download_file(
        config["urls"]["uniprot_fasta"],
        get_raw_path(config, "reference_fasta"),
        force=args.force_download,
    )
    structure_path = download_file(
        config["urls"]["alphafold_pdb"],
        get_raw_path(config, "alphafold_pdb"),
        force=args.force_download,
    )

    reference_sequence = read_fasta_sequence(fasta_path)
    clinvar_df = load_variant_summary(clinvar_path)
    cleaned = clean_clinvar_table(
        clinvar_df=clinvar_df,
        gene=config["gene"],
        reference_sequence=reference_sequence,
        domain_start=config["structure"]["domain_start"],
        domain_end=config["structure"]["domain_end"],
    )

    residue_df = extract_residue_table(structure_path, config["structure"]["chain_id"])
    residue_df.to_csv(get_processed_path(config, "cleaned_dataset").with_name("tp53_residue_table.csv"), index=False)

    cleaned = cleaned.merge(
        residue_df[["position", "structure_residue", "plddt", "x", "y", "z"]],
        on="position",
        how="inner",
    )
    cleaned = cleaned.loc[cleaned["structure_residue"] == cleaned["wild_type"]].copy()
    cleaned = assign_spatial_clusters(cleaned, residue_df)
    cleaned.to_csv(get_processed_path(config, "cleaned_dataset"), index=False)

    if args.foldseek:
        saprot_sequences = build_saprot_sequences(
            foldseek_path=args.foldseek,
            structure_path=structure_path,
            chain_id=config["structure"]["chain_id"],
            plddt_threshold=float(config["structure"]["plddt_mask_threshold"]),
        )
        save_saprot_sequences(saprot_sequences, get_processed_path(config, "saprot_sequences"))

    print(f"Saved cleaned dataset to {get_processed_path(config, 'cleaned_dataset')}")
    print(f"Number of variants: {len(cleaned)}")
    print(cleaned['label'].value_counts().to_dict())


if __name__ == "__main__":
    main()
