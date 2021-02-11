from Bio import SeqIO
import edlib
from sklearn.cluster import AffinityPropagation, DBSCAN, OPTICS
import time
import sys
import argparse
import pickle


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType(),
                        required=True)
    parser.add_argument('-o', '--output')
    parser.add_argument('-c', '--make_cluster_files',
                        dest='make_cluster_files', action='store_true')
    parser.add_argument('-d', '--file_dists', type=argparse.FileType("rb"))
    parser.set_defaults(make_cluster_files=False)
    return parser


def get_distance(str1, str2):
    return edlib.align(str1, str2)["editDistance"]


def get_matrix_of_dists(seq_records):
    matrix = [[0] * len(seq_records) for i in range(len(seq_records))]
    for i in range(len(seq_records)):
        for j in range(len(seq_records)):
            matrix[i][j] = get_distance(seq_records[i].seq, seq_records[j].seq)
    return matrix


if __name__ == '__main__':
    start_time = time.time()

    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])

    seqs = []
    for index, record in enumerate(SeqIO.parse(namespace.input, "fasta")):
        seqs.append(record)

    if namespace.file_dists:
        dists_matrix = pickle.load(namespace.file_dists)
    else:
        dists_matrix = get_matrix_of_dists(seqs)
        with open("dists " + namespace.input.name, 'wb') as fp:
            pickle.dump(dists_matrix, fp)

    print(time.time() - start_time, len(dists_matrix))
    # clustering = AgglomerativeClustering(affinity="precomputed",
    #                                      linkage="average").fit(
    #     dists_matrix)
    # clustering = AffinityPropagation(affinity="precomputed",
    #                                  random_state=0, max_iter=500).fit(dists_matrix)


    # clustering = DBSCAN(eps=5, metric="precomputed").fit(dists_matrix)
    clustering = OPTICS(min_samples=15, metric="precomputed").fit(dists_matrix)

    if namespace.make_cluster_files:
        files = []
        for i in set(clustering.labels_):
            filename = "cluster_" + str(i) + ".fa"
            file = open(filename, "w")
            file.close()

    if namespace.output is None:
        output_filename = "labels_" + namespace.input.name
    else:
        output_filename = namespace.output
    file_w = open(output_filename, "w")
    for i in range(len(seqs)):
        print(seqs[i].id, clustering.labels_[i], file=file_w)
        if namespace.make_cluster_files:
            filename = "cluster_" + str(clustering.labels_[i]) + ".fa"
            with open(filename, "a") as output_handle:
                SeqIO.write(seqs[i], output_handle, "fasta")

    file_w.close()
    print("Found", len(set(clustering.labels_)), "clusters:",
          *(sorted(list(set(clustering.labels_)))))
    print("Time:", time.time() - start_time)
