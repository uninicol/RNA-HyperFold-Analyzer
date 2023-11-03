import warnings

from Bio import SeqIO

from incidence_producers.vienna_incidence_producer import ViennaIncidenceProducer
from rna_stats.forna_rna_stats import FornaRnaStats
from temp_classes.rna_folder import RNAFolder

# def timefunc(function, *args):
#     start_time = time.perf_counter()
#     function(*args)
#     print(f"--- {time.perf_counter() - start_time:.2f} seconds ---")


# def temperature_dotbracket(parallel=False):
#     bracket_set = set()
#     if parallel:
#         with Pool() as pool:
#             for result in pool.imap(
#                 folder.temperature_folding, range(start_temp, end_temp + 1)
#             ):
#                 bracket_set.add(result)
#     else:
#         for temperature in range(start_temp, end_temp + 1):
#             bracket_set.add(folder.temperature_folding(temperature))
#     print(f"numero di folding rilevati:\t{len(bracket_set)}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")  # non vengono scritti i warning di pandas

    fasta_file = "examples/small_sequence.fasta"
    sequence = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
    # thg = TemperatureHypergraph(ViennaIncidenceProducer(RNAFolder(sequence)))
    # thg.insert_temperature(5)
    # thg.insert_temperature(100)
    # fasta_file = "examples/sequence (2).fasta"
    # start_temp = 94
    # end_temp = 98
    # start_temp = max(start_temp, -273)
    # sequence = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
    # folder = RNAFolder(sequence)

    # print(f"{start_temp=}, {end_temp=}")
    # print(f"lunghezza sequenza:\t{len(sequence)}")
    # print("parallelo:")
    # timefunc(temperature_dotbracket, True)
    # print("sequenziale:")
    # timefunc(temperature_dotbracket, False)

    for producer in [
        # FornaIncidenceProducer(forna_file),
        # ForgiIncidenceProducer(forna_file),
        ViennaIncidenceProducer(RNAFolder(sequence)),
    ]:
        rna_stats = FornaRnaStats(producer, 0)
        rna_stats.plot_hypergraph()
        rna_stats = FornaRnaStats(producer, 100)
        rna_stats.plot_hypergraph()
    #     print(rna_stats.secondary_structures())
    #     print(f"modularity={rna_stats.modularity()}")
    #     print("-" * 70)
    # #     #
    #     rna_stats.plot_partitions_conductance()
    # #     # # for partition, conductance in enumerate(rna_stats.partitions_conductance()):
    # #     # #     print(f"{partition=}, {conductance=}")
    # #     # # print("-" * 70)
    # #     # # for edge, centrality in rna_stats.n_between_centrality(n=1).items():
    # #     # #     print(f"{edge=}, {centrality=}")
    #     rna_stats.plot_n_between_centrality(n=1)
    #     print(rna_stats.partitions())
    #     print(rna_stats.secondary_structures())
    #     # rna_stats.plot_n_between_centrality(n=2)
pass
