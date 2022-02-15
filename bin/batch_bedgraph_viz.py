#!/usr/bin/env python3

import boto3
import pandas as pd
import os
import sys
import logging


def get_file_names(s3path, suffix="txt"):
    """Get a list of s3paths given certain restrictions on prefix and suffix

    Args:
        s3path (str): Only fetch keys that start with this prefix (folder name).
        suffix (str, optional): Only fetch keys that end with this suffix (extension). Defaults to "txt".

    Returns:
        list: all the file names in an S3 bucket folder.
    """
    s3_client = boto3.client("s3")
    bucket_name, prefix = split_s3_path(s3path)
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
    objs = response["Contents"]

    while response["IsTruncated"]:
        response = s3_client.list_objects_v2(
            Bucket=bucket_name,
            Prefix=prefix,
            ContinuationToken=response["NextContinuationToken"],
        )
        objs.extend(response["Contents"])

    logging.info(f"Sifting through {len(objs)} files ... ")

    shortlisted_files = list()
    if suffix == "":
        shortlisted_files = [obj["Key"] for obj in objs]
        total_size_bytes = sum([obj["Size"] for obj in objs])
    else:
        logging.info(f"Limiting search by suffix: {suffix}")
        shortlisted_files = [obj["Key"] for obj in objs if obj["Key"].endswith(suffix)]
        total_size_bytes = sum(
            [obj["Size"] for obj in objs if obj["Key"].endswith(suffix)]
        )

    logging.info(
        f"Found {len(shortlisted_files)} files, totalling about {total_size_bytes/1e9:,.3f} Gb."
    )
    return shortlisted_files


def split_s3_path(s3_path):

    try:
        split_path = s3_path.replace("s3://", "").split("/")
    except:
        logging.error(s3_path)
        sys.exit(1)

    bucket = split_path[0]
    key = "/".join(split_path[1:])
    return (bucket, key)


def download_from_s3(bucket_name, object_keys, local_dir):
    downloaded_paths = list()
    s3 = boto3.client("s3")
    for object_key in object_keys:
        local_file_path = os.path.join(local_dir, os.path.basename(object_key))
        if not os.path.exists(local_file_path):
            s3.download_file(bucket_name, object_key, local_file_path)
        downloaded_paths.append(local_file_path)
    return downloaded_paths


def download_bedgraphs(df, suffix, local_bedgraph_dir, bucket_name):
    bedgraph_files_series = df["output_paths"].apply(
        lambda x: get_file_names(x, suffix)
    )

    local_bedgraph_files = list()
    for bedgraph_file in bedgraph_files_series:
        local_bedgraph_files.extend(
            download_from_s3(bucket_name, bedgraph_file, local_bedgraph_dir)
        )

    return local_bedgraph_files


def create_igv_stub_per_genome(
    genome_name, genome_fasta, genome_bedgraph_files, genome_images_path, suffix
):
    ## for each genome
    stub = list()
    # buffer for next genome
    stub.append("##")
    # new
    stub.append("new")
    # genome /Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/fasta/Acidaminococcus-sp-D21.fna
    stub.append(f"genome {genome_fasta}")
    # load  /Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/bedgraphs/Dropout_Acidaminococcus-sp-D21_t24_MM_replicate1--Acidaminococcus-sp-D21.primary.bedgraph
    # load  /Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/bedgraphs/Dropout_Acidaminococcus-sp-D21_t48_MM_replicate1--Acidaminococcus-sp-D21.primary.bedgraph
    # load  /Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/bedgraphs/Dropout_Acidaminococcus-sp-D21_t48_SACC_replicate1--Acidaminococcus-sp-D21.primary.bedgraph
    for bg_file in genome_bedgraph_files:
        stub.append(f"load {bg_file}")
    # snapshotDirectory /Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/IGV_snapshots/Acidaminococcus-sp-D21/
    stub.append(f"snapshotDirectory {genome_images_path}")
    # snapshot Acidaminococcus-sp-D21.svg
    stub.append(f"snapshot {genome_name}.{suffix}.svg")
    return stub


if __name__ == "__main__":
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    metadata_file = "/Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/StrainDropouts.metadata.csv"
    parent_s3path = "s3://genomics-workflow-core/Results/ninjamap2/20220201_Strain_Dropout_Verification"
    file_suffix = "primary.bedgraph"
    file_suffix_clean = file_suffix.replace(".", "_")
    output_dir = "/Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/StrainDropouts"
    fasta_dir = (
        "/Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/fasta"
    )
    batch_script_filename = "/Users/sunit.jain/Research/Alice/revisions/NM_Verify/SCv1_2_20211008/StrainDropouts/20220202_StrainDropouts.igv_batch_script.txt"
    os.makedirs(output_dir, exist_ok=True)

    parent_bucket_name, _ = split_s3_path(parent_s3path)

    metadata_df = pd.read_csv(metadata_file)
    metadata_df["output_paths"] = metadata_df["sample_name"].apply(
        lambda x: f"{parent_s3path}/{x}"
    )

    grouped_metadata_df = metadata_df.groupby(["source_genome"])
    batch_script_cmds = list()
    for name, group_df in grouped_metadata_df:
        logging.info(f"Processing {name}")
        genome_fasta = os.path.join(fasta_dir, f"{name}.fna")

        genome_parent_path = os.path.join(output_dir, name)
        genome_bedgraphs_path = os.path.join(genome_parent_path, "bedgraphs")
        genome_images_path = os.path.join(output_dir, "images", file_suffix_clean)

        os.makedirs(genome_bedgraphs_path, exist_ok=True)
        os.makedirs(genome_images_path, exist_ok=True)

        genome_bedgraph_files = download_bedgraphs(
            group_df, file_suffix, genome_bedgraphs_path, parent_bucket_name
        )

        if len(genome_bedgraph_files) > 0:
            batch_script_cmds.extend(
                create_igv_stub_per_genome(
                    name,
                    genome_fasta,
                    genome_bedgraph_files,
                    genome_images_path,
                    file_suffix_clean,
                )
            )

    cmds = "\n".join(batch_script_cmds)
    # logging.info(cmds)
    with open(batch_script_filename, "w") as batch_script:
        batch_script.write(cmds)
