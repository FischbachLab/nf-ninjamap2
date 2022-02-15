#!/usr/bin/env python3

import json
import logging
import os

# import pickle
import sys
from time import perf_counter as timer

import numpy as np
import pysam

# from matplotlib import pyplot as plt


def human_time(time):
    time = abs(time)
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    time_str = format("%02d:%02d:%02d:%02d" % (day, hour, minutes, seconds))
    return time_str


def is_perfect_alignment(aln):
    edit_dist = dict(aln.tags)["NM"]
    query_len = aln.query_length
    ref_start = aln.reference_start
    ref_end = aln.reference_end

    # https://www.biostars.org/p/106126/
    return (edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end))


# Each list of alignments looks like this
# Forward read - Acidaminococcus-sp-D21-fna_Node_75-4258 :
# fwd_list = [('Clostridium-scindens-ATCC-35704-fna_Node_7', 1272, 1422), ('Anaerotruncus-colihominis-DSM-17241-fna_Node_9', 376197, 376347), ('Dorea-formicigenerans-ATCC-27755-fna_Node_2', 334213, 334363), ('Anaerostipes-caccae-DSM-14662-fna_Node_9', 177840, 177990), ('Ruminococcus-gnavus-ATCC-29149-MAF-NJ26-fna_Node_0', 1385916, 1386066), ('Coprococcus-eutactus-ATCC-27759-fna_Node_3', 44782, 44932), ('Clostridium-nexile-DSM-1787-fna_Node_6', 78316, 78466), ('Acidaminococcus-sp-D21-fna_Node_75', 42945, 43095), ('Clostridium-sp-M62-1-fna_Node_24', 245079, 245229)]

# Reverse read - Acidaminococcus-sp-D21-fna_Node_75-4258 :
# rev_list = [('Clostridium-scindens-ATCC-35704-fna_Node_7', 1612, 1762), ('Anaerotruncus-colihominis-DSM-17241-fna_Node_9', 376537, 376687), ('Dorea-formicigenerans-ATCC-27755-fna_Node_2', 333873, 334023), ('Anaerostipes-caccae-DSM-14662-fna_Node_9', 178180, 178330), ('Ruminococcus-gnavus-ATCC-29149-MAF-NJ26-fna_Node_0', 1385576, 1385726), ('Coprococcus-eutactus-ATCC-27759-fna_Node_3', 44442, 44592), ('Clostridium-nexile-DSM-1787-fna_Node_6', 77976, 78126), ('Acidaminococcus-sp-D21-fna_Node_75', 42605, 42755), ('Clostridium-sp-M62-1-fna_Node_24', 245419, 245569)]
def get_common_genomes(fwd_list, rev_list, bin_map):
    fwd_genomes = {bin_map[contig] for contig, _, _ in fwd_list}
    rev_genomes = {bin_map[contig] for contig, _, _ in rev_list}

    common_genomes = fwd_genomes.intersection(rev_genomes)

    common_aln_fwd = [
        (contig, start, stop)
        for contig, start, stop in fwd_list
        if bin_map[contig] in common_genomes
    ]
    common_aln_rev = [
        (contig, start, stop)
        for contig, start, stop in rev_list
        if bin_map[contig] in common_genomes
    ]

    common_aln = common_aln_fwd + common_aln_rev
    return common_genomes, common_aln


def add_to_bin(read_name, vote_bin_aln, vote_bin_reads, contig_tuple, bin_map):
    contig = contig_tuple[0]
    genome = bin_map[contig]
    if genome in vote_bin_aln:
        vote_bin_aln[genome].append(contig_tuple)
    else:
        vote_bin_aln.update({genome: [contig_tuple]})

    if read_name in vote_bin_reads:
        vote_bin_reads[read_name].add(genome)
    else:
        vote_bin_reads.update({read_name: {genome}})
    return


## primary_reads and escrow_reads look like this: read_name --> set of genomes
# {'Acidaminococcus-sp-D21-fna_Node_75-4258': {'Acidaminococcus-sp-D21-fna',
#   'Anaerostipes-caccae-DSM-14662-fna',
#   'Anaerotruncus-colihominis-DSM-17241-fna',
#   'Clostridium-nexile-DSM-1787-fna',
#   'Clostridium-scindens-ATCC-35704-fna',
#   'Clostridium-sp-M62-1-fna',
#   'Coprococcus-eutactus-ATCC-27759-fna',
#   'Dorea-formicigenerans-ATCC-27755-fna',
#   'Ruminococcus-gnavus-ATCC-29149-MAF-NJ26-fna'},
#  'Acidaminococcus-sp-D21-fna_Node_75-5856': { ... },
#     ...}
def aggregate_votes_from_reads(bin_reads, pair_contribution, library_size, priors=None):
    votes = dict()
    read_votes = dict()

    default_weight = 0
    num_voted_reads = 0

    if priors is None:
        priors = dict()
        default_weight = 1

    for read_name in bin_reads:
        list_of_aln_genomes = list()
        read_vote_contribution = 0
        sum_genome_weights = 0
        expected_total_vote = (
            pair_contribution[read_name] if read_name in pair_contribution else 1
        )
        read_votes.update({read_name: {}})

        for genome in bin_reads[read_name]:
            list_of_aln_genomes.append(genome)
            genome_weight = priors[genome] if genome in priors else default_weight
            sum_genome_weights += genome_weight

        if sum_genome_weights == 0:
            # might happen if none of the genomes had a prior.
            logging.debug(
                f"Skipping {read_name} with a possible total vote contribution of {expected_total_vote} because it aligned to the follwoing genomes, none of which had a prior: {list_of_aln_genomes}"
            )
            continue

        num_voted_reads += 1

        for genome in bin_reads[read_name]:
            genome_weight = priors[genome] if genome in priors else default_weight
            adjusted_vote = expected_total_vote * (genome_weight / sum_genome_weights)
            read_votes[read_name].update({genome: adjusted_vote})
            read_vote_contribution += adjusted_vote
            if genome in votes:
                votes[genome] += adjusted_vote
            else:
                votes.update({genome: adjusted_vote})

        assert (
            int(round(read_vote_contribution)) == expected_total_vote
        ), logging.critical(
            f"Voter fraud! {read_name} voted {read_vote_contribution} times! Should have been {expected_total_vote}."
        )

    # normalize by library size
    for genome, num_votes in votes.items():
        votes[genome] = num_votes / library_size

    return num_voted_reads, votes, read_votes


# vote_bin looks like this:
# {'Acidaminococcus-sp-D21-fna':
#  [('Acidaminococcus-sp-D21-fna_Node_0', 2, 152),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 355, 505),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 12, 162),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 371, 521),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 14, 164),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 366, 516),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 18, 168),
# ...]
# }
# genome_size = {'Acidaminococcus-sp-D21-fna': {'Acidaminococcus-sp-D21-fna_Node_0':2305309 }}
def get_genome_coverage_array(vote_bin, bin_stats, vote_value=1):
    genome_coverage = dict()
    for genome in vote_bin:
        genome_contigs_size = bin_stats[genome]
        genome_cov_dict = dict()

        for contig, start, end in vote_bin[genome]:
            if contig not in genome_cov_dict:
                contig_arr = np.zeros((genome_contigs_size[contig] + 1))
                genome_cov_dict.update({contig: contig_arr})

            contig_arr = genome_cov_dict[contig]
            contig_arr[start:end] += vote_value

        # Better to make just one adjustment rather than updating each entry.
        [np.delete(contig_cov_arr, 0) for _, contig_cov_arr in genome_cov_dict.items()]
        genome_coverage.update({genome: genome_cov_dict})
    return genome_coverage


## Genome coverage dictionary looks like this for each genome:
# {'Clostridium-scindens-ATCC-35704-fna_Node_7': array([0., 0., 0., ..., 0., 0., 0.]),
#   'Clostridium-scindens-ATCC-35704-fna_Node_11': array([0., 0., 0., ..., 0., 0., 0.]),
#   'Clostridium-scindens-ATCC-35704-fna_Node_40': array([0., 0., 0., ..., 0., 0., 0.])},

#  {'Anaerotruncus-colihominis-DSM-17241-fna_Node_9': array([0., 0., 0., ..., 0., 0., 0.]),
#   'Anaerotruncus-colihominis-DSM-17241-fna_Node_17': array([0., 0., 0., ..., 0., 0., 0.])},
def genome_cov_stats(genome, genome_dict, bin_stats):
    genome_stats = dict()
    genome_stats.update({genome: {}})
    for contig, contig_cov_arr in genome_dict.items():
        # Calculate stats
        genome_stats[genome].update(
            {
                contig: {
                    "Length": bin_stats[genome][contig],
                    "Mean_Coverage": np.mean(contig_cov_arr),
                    "Stdev_Coverage": np.std(contig_cov_arr),
                    "Covered_bases": np.sum(contig_cov_arr > 0),
                }
            }
        )
    return genome_stats


def read_binmap_file(binmap_file):
    genome_size = dict()
    contig_to_genome_map = dict()
    genome_idx = 0
    contig_idx = 6
    length_idx = 7
    with open(binmap_file, "r") as binmap:
        for idx, line in enumerate(binmap):
            if idx == 0:
                continue
            content = line.rstrip().split(",")
            # genome = content[genome_idx].replace("-fna", "")
            # contig = content[contig_idx].replace("-fna", "")
            genome = content[genome_idx]
            contig = content[contig_idx]
            length = int(content[length_idx])

            if genome not in genome_size:
                #             genome_size.update({genome:{"Length":{},"Mean_Coverage":{},"Stdev_Coverage":{}}})
                genome_size.update({genome: {}})

            genome_size[genome].update({contig: length})

            if contig not in contig_to_genome_map:
                contig_to_genome_map.update({contig: genome})
    # genome_size["Acidaminococcus-sp-D21-fna"]["Acidaminococcus-sp-D21-fna_Node_0"]
    return genome_size, contig_to_genome_map


def get_reads_with_perfect_alignments(local_bam, debug=False):
    fwd_reads = dict()
    rev_reads = dict()
    bamfile = pysam.AlignmentFile(local_bam, mode="rb")
    for aln in bamfile.fetch(until_eof=True):
        # if percent id and alignment length is not 'perfect', move on
        if not is_perfect_alignment(aln):
            continue

        read_name = str(aln.qname)
        if debug:
            read_source = read_name.split("_")[0]
        contig_name = str(aln.reference_name)

        # num_aln += 1

        start, stop = sorted((int(aln.reference_start), int(aln.reference_end)))

        # if read is FWD/REV:
        # add it to the appropriate fwd_reads/rev_reads dict with the name of
        # the contig that it aligned to, and the sorted start and end position
        # of the alignment.
        if aln.is_read1:
            if read_name in fwd_reads:
                fwd_reads[read_name].append((contig_name, start, stop))
            else:
                fwd_reads.update({read_name: [(contig_name, start, stop)]})
        else:
            if read_name in rev_reads:
                rev_reads[read_name].append((contig_name, start, stop))
            else:
                rev_reads.update({read_name: [(contig_name, start, stop)]})

    return (fwd_reads, rev_reads)


def read_partitioning(fwd_reads, rev_reads, contig_map):
    primary_aln = dict()
    escrow_aln = dict()
    pair_contribution = dict()

    primary_reads = dict()
    escrow_reads = dict()

    # for each fwd read
    for read_name in fwd_reads:
        # get a list of alignments for fwd reads;
        # each alignment is a tuple with (contig name, start position of aln on contig, end position of aln on contig)
        list_of_fwd_aln = fwd_reads[read_name]

        # if mate does not have a perfect aln, put this fwd read in escrow
        if read_name not in rev_reads:
            [
                add_to_bin(
                    read_name, escrow_aln, escrow_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_fwd_aln
            ]
            continue

        # get a list of alignments for fwd reads;
        list_of_rev_aln = rev_reads[read_name]

        # Since the pair had perfect alignment,
        # process it here and avoid reprocessing it later
        # by deleting it from the rev_reads aln dictionary
        del rev_reads[read_name]

        # identify if there are any common genomes that the pair aligns to
        # and aggregate all alignments.
        # If the pair align to different genomes, identify the common genomes and
        # keep only the alignments to the common genomes and discard the remaining
        # If there are no common genomes, ignore all alignment, they will all be added to the
        # escrow bin later
        (set_of_common_genomes, list_of_common_aln,) = get_common_genomes(
            list_of_fwd_aln, list_of_rev_aln, contig_map
        )

        # if the pair has exactly one genome in common,
        # this is the BEST case and use this to calculate priors later
        # for now, add both reads to the primary bin
        if len(set_of_common_genomes) == 1:
            # add_to_primary_bin:
            #     - assign 2 votes to the genome
            #     - record the votes
            [
                add_to_bin(
                    read_name, primary_aln, primary_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_common_aln
            ]
            pair_contribution.update({read_name: 2})
        # if the pair has more than one genome in common,
        # add these reads and just the common genome alignments to the escrow bin
        elif len(set_of_common_genomes) > 1:
            [
                add_to_bin(
                    read_name, escrow_aln, escrow_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_common_aln
            ]
            pair_contribution.update({read_name: 2})
        # if the pair has no genome in common,
        # add these reads and all their alignments to the escrow bin
        elif len(set_of_common_genomes) == 0:
            [
                add_to_bin(
                    read_name, escrow_aln, escrow_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_fwd_aln
            ]
            [
                add_to_bin(
                    read_name, escrow_aln, escrow_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_rev_aln
            ]

    # Since all rev reads that had a match in the fwd
    # have already been removed, just in case there are extra
    # put them into escrow with all their alignments.
    if len(rev_reads) > 0:
        for read_name in rev_reads:
            list_of_rev_aln = rev_reads[read_name]
            [
                add_to_bin(
                    read_name, escrow_aln, escrow_reads, contig_tuple, contig_map
                )
                for contig_tuple in list_of_rev_aln
            ]

    return (primary_aln, escrow_aln, primary_reads, escrow_reads, pair_contribution)


def tally_votes(primary, escrow, pair_contribution, library_size):
    tally = dict()
    (
        num_voted_primary,
        primary_votes,
        primary_votes_by_reads,
    ) = aggregate_votes_from_reads(primary, pair_contribution, library_size)

    # use the primary votes as priors
    num_voted_escrow, escrow_votes, escrow_votes_by_reads = aggregate_votes_from_reads(
        escrow, pair_contribution, library_size=library_size, priors=primary_votes,
    )

    for genome, primary_vote in primary_votes.items():
        escrow_vote = escrow_votes[genome] if genome in escrow_votes else 0
        tally.update({genome: (primary_vote + escrow_vote)})

    return (
        num_voted_primary,
        num_voted_escrow,
        primary_votes_by_reads,
        escrow_votes_by_reads,
        tally,
    )


def save_dict_as_json(dictionary, filename):
    # create json object from dictionary
    content = json.dumps(dictionary)

    # open file for writing, "w"
    f = open(filename, "w")
    # write json object to file
    f.write(content)
    # close file
    f.close()


# def save_dict_as_pickle(dictionary, filename):
#     f = open(filename, "wb")
#     # write the python object (dict) to pickle file
#     pickle.dump(dictionary, f)

#     # close file
#     f.close()


# primary_aln and escrow_aln looks like this:
# {'Acidaminococcus-sp-D21-fna':
#  [('Acidaminococcus-sp-D21-fna_Node_0', 2, 152),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 355, 505),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 12, 162),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 371, 521),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 14, 164),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 366, 516),
#   ('Acidaminococcus-sp-D21-fna_Node_0', 18, 168),
# ...]
# }

## Genome coverage dictionary looks like this for each genome:
# {'Clostridium-scindens-ATCC-35704-fna_Node_7': array([0., 0., 0., ..., 0., 0., 0.]),
#   'Clostridium-scindens-ATCC-35704-fna_Node_11': array([0., 0., 0., ..., 0., 0., 0.]),
#   'Clostridium-scindens-ATCC-35704-fna_Node_40': array([0., 0., 0., ..., 0., 0., 0.])},
def save_aln_as_bedgraph(read_aln, coverage, target_genome, prefix, read_type):
    #     https://genome.ucsc.edu/goldenPath/help/bedgraph.html
    # one file for each genome
    with open(f"{prefix}-{target_genome}.{read_type}.bedGraph", "w") as bg:
        header = f'track type=bedGraph name="{prefix}" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'
        bg.write(f"{header}\n")
        for contig, start, stop in sorted(read_aln[target_genome]):
            mean_coverage = np.mean(coverage[target_genome][contig])
            bg.write(f"{contig}\t{start}\t{stop}\t{mean_coverage}\n")
    return


## Coverage array dict looks like this
# {'Clostridium-scindens-ATCC-35704-fna': {
#                  'Clostridium-scindens-ATCC-35704-fna_Node_7': np.array([0., 0., 1., 1., 1., 0.]),
#                  'Clostridium-scindens-ATCC-35704-fna_Node_11': np.array([2., 2., 1., 1., 0., 0.]),
#                  'Clostridium-scindens-ATCC-35704-fna_Node_40': np.array([3., 0., 0., 0., 0., 0.])
# },
#              'Anaerotruncus-colihominis-DSM-17241-fna': {
#                  'Anaerotruncus-colihominis-DSM-17241-fna_Node_9': np.array([1., 0., 0., 0., 0., 0.]),
#                  'Anaerotruncus-colihominis-DSM-17241-fna_Node_17': np.array([2., 0., 0., 0., 0., 0.])
#              }
#             }
def save_cov_array_as_bedgraph(coverage, target_genome, prefix, read_type, plot=False):
    if target_genome not in coverage:
        logging.info(
            f"Genome '{target_genome}' does not appear to have any coverage of type '{read_type}'. Skipping ..."
        )
        return

    filename = f"{prefix}--{target_genome}"
    with open(f"{filename}.{read_type}.bedgraph", "w") as bg:
        header = f'track type=bedGraph name="{prefix}" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'
        bg.write(f"{header}\n")
        for contig in coverage[target_genome]:
            prev_base_cov = -1
            cov_change_start = -1
            contig_len = len(coverage[target_genome][contig])
            # if plot:
            #     plot_coverage(
            #         coverage[target_genome][contig], contig, f"{filename}--{contig}"
            #     )
            for idx, current_coverage in enumerate(coverage[target_genome][contig]):
                if prev_base_cov == -1:
                    prev_base_cov = current_coverage
                    cov_change_start = 0
                    continue

                if current_coverage != prev_base_cov:
                    changed_on = idx - 1
                    bg.write(
                        f"{contig}\t{cov_change_start}\t{changed_on}\t{prev_base_cov}\n"
                    )
                    prev_base_cov = current_coverage
                    cov_change_start = idx

                if idx == contig_len - 1:
                    bg.write(
                        f"{contig}\t{cov_change_start}\t{idx}\t{current_coverage}\n"
                    )
        return


# def plot_coverage(data, contig, title, prefix):
#     plt.title(f"{title}")
#     plt.ylabel(f"Coverage")
#     plt.xlabel(f"{contig}")
#     plt.plot(data)
#     plt.savefig(f"{prefix}.png", dpi=100, bbox_inches="tight")
#     plt.close()
#     return


##############################################################################
# Main
##############################################################################

if __name__ == "__main__":
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.DEBUG,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    start = timer()
    logging.info("Started")
    # bamfile_name = "s3://genomics-workflow-core/Results/NinjaIndex/DaisyTest4/db/bowtie2_mapping/tmp_20210827_213525/Sync/bowtie2/Collinsella-intestinalis-SH0000052-MAF-2.R1.fastq.gz_vs_db_all_genomes.name_sorted.markdup.bam"
    # bamfile_idx_name = "s3://genomics-workflow-core/Results/NinjaIndex/DaisyTest4/db/bowtie2_mapping/tmp_20210827_213525/Sync/bowtie2/Collinsella-intestinalis-SH0000052-MAF-2.R1.fastq.gz_vs_db_all_genomes.name_sorted.markdup.bam.bai"
    # local_path = "/Users/sunit.jain/Research/NinjaMap/test_data/ninjaindex_bams"

    # local_bam = "/Users/sunit.jain/Research/NinjaMap/Database/SCv1_2_202110/Bams/Acidaminococcus-sp-D21-fna.R1.fastq.gz_vs_db_all_genomes.name_sorted.markdup.bam"
    # binmap_file = "/Users/sunit.jain/Research/NinjaMap/Database/SCv1_2_202110/SCv1_2.ninjaIndex.binmap.csv"
    # library_size = 0
    # prefix = "Acidaminococcus-sp-D21-fna"

    local_bam = sys.argv[1]
    binmap_file = sys.argv[2]
    library_size = int(sys.argv[3])
    prefix = sys.argv[4]
    target_genome = sys.argv[5]

    # Output
    primary_reads_list = f"{prefix}.primary_reads.tsv"
    primary_alignments_file = f"{prefix}.primary_aln.json"
    escrow_alignments_file = f"{prefix}.escrow_aln.json"
    escrow_coverage_file = f"{prefix}.escrow_coverage.json"
    escrow_coverage_stats_file = f"{prefix}.escrow_coverage_stats.json"
    primary_coverage_file = f"{prefix}.primary_coverage.json"
    primary_coverage_stats_file = f"{prefix}.primary_coverage_stats.json"
    read_stats_file = f"{prefix}.reads_stats.json"
    read_fractions_file = f"{prefix}.read_fractions.json"
    primary_votes_by_reads_file = f"{prefix}.primary_votes_by_reads.json"
    escrow_votes_by_reads_file = f"{prefix}.escrow_votes_by_reads.json"

    read_stats = dict()
    read_stats.update({"num_individual_reads_input": library_size})

    # Get some basic information about the database
    genome_size_dict, contig_map = read_binmap_file(binmap_file)
    logging.info(
        f"Read binmap file containing {len(genome_size_dict)} genomes with {len(contig_map)} contigs ..."
    )

    # filter for perfect individual read alignments
    # Parse the bam file one alignment at a time
    # Just looking for reads with perfect alignments, not
    # interested in whether the mate aligned or not, just yet.
    logging.info(f"Reading BAM file to identify reads with 'perfect' alignments ...")
    fwd_reads, rev_reads = get_reads_with_perfect_alignments(local_bam, debug=False)
    num_reads_with_perfect_aln = len(fwd_reads) + len(rev_reads)
    read_stats.update({"num_reads_with_perfect_aln": num_reads_with_perfect_aln})

    # TODO @sunitj
    # This is a place holder for testing.
    # Should be replaced with the library size parameter,
    # so the user/previous step in the pipeline provides this information.
    lib_size = library_size if library_size > 0 else num_reads_with_perfect_aln

    logging.info(
        f"Found {num_reads_with_perfect_aln} reads (about {num_reads_with_perfect_aln/2} pairs) with 'perfect' alignments."
    )

    # filter for read pairs with perfect alignments
    # Of the reads that were perfectly aligned, how many pairs were aligned to
    # the same genome? Set them aside as Primary and everything else as Escrow.
    (
        primary_aln,
        escrow_aln,
        primary_reads,
        escrow_reads,
        pair_contribution,
    ) = read_partitioning(fwd_reads, rev_reads, contig_map)
    del fwd_reads
    del rev_reads

    num_primary_reads = len(primary_reads)
    read_stats.update({"num_reads_with_primary_aln": num_primary_reads})
    logging.info(
        f"Found {num_primary_reads} read pairs with 'perfect' alignments to a single genome."
    )

    num_escrow_reads = len(escrow_reads)
    read_stats.update({"num_reads_with_escrow_aln": num_escrow_reads})
    logging.info(
        f"Found {num_escrow_reads} read pairs with 'perfect' alignments to more than one genome."
    )

    logging.info(
        f"Overall {len(pair_contribution)} read pairs agreed on at least one genome"
    )

    logging.info(f"Saving alignment JSONs ... ")
    save_dict_as_json(primary_aln, primary_alignments_file)
    save_dict_as_json(escrow_aln, escrow_alignments_file)
    logging.info(f" Done.")

    logging.info(f"Tallying votes for all genomes in the database ...")
    # Aggregate and normalize all the votes.
    (
        num_voted_primary,
        num_voted_escrow,
        primary_votes_by_reads,
        escrow_votes_by_reads,
        read_fractions,
    ) = tally_votes(primary_reads, escrow_reads, pair_contribution, lib_size)
    del primary_reads
    del escrow_reads
    logging.info(f" Done.")
    read_stats.update({"num_reads_with_primary_aln_voted": num_voted_primary})
    logging.info(
        f"Out of {num_primary_reads}, {num_voted_primary} voted as primary reads"
    )
    read_stats.update({"num_reads_with_escrow_aln_voted": num_voted_escrow})
    logging.info(f"Out of {num_escrow_reads}, {num_voted_escrow} voted as escrow reads")

    logging.info(f"Saving voting JSONs ... ")
    save_dict_as_json(read_fractions, read_fractions_file)
    save_dict_as_json(primary_votes_by_reads, primary_votes_by_reads_file)
    save_dict_as_json(escrow_votes_by_reads, escrow_votes_by_reads_file)
    logging.info(f" Done.")

    logging.info(f"Getting genome coverage from Primary reads ... ")
    primary_coverage_arr = get_genome_coverage_array(primary_aln, genome_size_dict, 1)
    primary_coverage_stats_list = [
        genome_cov_stats(genome, genome_dict, genome_size_dict)
        for genome, genome_dict in primary_coverage_arr.items()
    ]
    logging.info(f" Done.")

    # logging.info(f"Saving primary coverage files ... ")
    # save_dict_as_pickle(primary_coverage_arr, primary_coverage_file)
    # save_dict_as_pickle(primary_coverage_stats_list, primary_coverage_stats_file)
    # logging.info(f" Done.")

    logging.info(f"Getting genome coverage from Escrow reads ... ")
    escrow_coverage_arr = get_genome_coverage_array(escrow_aln, genome_size_dict, 1)
    escrow_coverage_stats_list = [
        genome_cov_stats(genome, genome_dict, genome_size_dict)
        for genome, genome_dict in escrow_coverage_arr.items()
    ]
    logging.info(f" Done.")

    logging.info(f"Saving {target_genome} genome coverage to bedGraph format ...")
    # save_aln_as_bedgraph(
    #     primary_aln, primary_coverage_arr, target_genome, prefix, read_type="primary"
    # )
    save_cov_array_as_bedgraph(
        primary_coverage_arr, target_genome, prefix, read_type="primary"
    )
    save_cov_array_as_bedgraph(
        escrow_coverage_arr, target_genome, prefix, read_type="escrow"
    )
    logging.info(f" Done.")

    # logging.info(f"Saving escrow coverage files ... ")
    # save_dict_as_pickle(escrow_coverage_arr, escrow_coverage_file)
    # save_dict_as_pickle(escrow_coverage_stats_list, escrow_coverage_stats_file)
    # logging.info(f" Done.")

    logging.info(f"Saving read stats ... ")
    save_dict_as_json(read_stats, read_stats_file)
    logging.info(f" Done.")
    end = timer()
    logging.info("Completed in %s (d:h:m:s)", human_time(end - start))
    logging.info(f"Pipeline finished without errors! Huzzah!")

