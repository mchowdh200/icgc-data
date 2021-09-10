import sys
import random
import numpy as np
from pprint import pprint
# from gensim.models import HdpModel, LdaModel, TfidfModel
# from gensim.corpora import Dictionary
import tomotopy as tp
import matplotlib.pyplot as plt
from wordcloud import WordCloud

def construct_corpus(input_files: list[str], feature_column: int,
                     count_column: int) -> list[list[str]]:
    """
    Contstruct a list of samples where each sample will contain
    a list of features present in that sample (stored as strings).
    Returns 
    """

    # first index is the "document level"
    # second will be the "word level"
    samples = []
    for bed in input_files:
        S = []
        with open(bed) as f:
            for line in f:
                if int((x:=line.rstrip().split())[count_column]) > 0:
                    S.extend([x[feature_column]] * int(x[count_column]))
            samples.append(S)
    return samples

def train_hdp(corpus, initial_k=10, term_weight=tp.TermWeight.PMI , iterations=1000):
    """
    Train a heirarchical dirichlet process topic model
    """

    hdp = tp.HDPModel(tw=term_weight,
                      gamma=1,
                      alpha=0.1,
                      initial_k=initial_k,
                      seed=1000)

    # add samples to model
    random.seed(1000)
    random.shuffle(corpus)
    for c in corpus:
        hdp.add_doc(c)

    # discard the first N samples
    hdp.burn_in = 1000
    hdp.train(0)

    for i in range(0, iterations+1, 100):
        hdp.train(100)
        print(f'{i = }\tlog-likelihood = {hdp.ll_per_word}\ttopics = {hdp.live_k}')
        
    return hdp


def get_topics(hdp) -> list[list[tuple[str, float]]]:
    # sort topic by importance (descending)
    sorted_topics = [i for i, _ in sorted(
        enumerate(hdp.get_count_by_topics()),
        key=lambda x: x[1], reverse=True)]

    # return only the live topics along with their word/probability arrays
    return [hdp.get_topic_words(i, top_n=15) for i in sorted_topics
            if hdp.is_live_topic(i)]


def closestDivisors(n):
    """
    dumb way I thought of to create a good grid of subplots
    """
    a = np.round(np.sqrt(n))
    while n%a > 0:
        a -= 1
    return a, n//a


def visualize_topics(topics):
    """
    Generate word cloud visualizations for each topic in
    the provided topic model (saved as png).
    """
    N = len(topics)
    nrows, ncols = closestDivisors(N)
    for i, topic in enumerate(topics):
        freqs = {word : freq for word, freq in topic}
        # wc = WordCloud(background_color="black", max_words=100)
        # wc.generate_from_frequencies(freqs)
        plt.subplot(nrows, ncols, i+1)
        plt.barh(range(len(freqs)), freqs.values())
        plt.yticks(range(len(freqs)), list(freqs.keys()))
        # plt.imshow(wc, interpolation="bilinear")
        # plt.axis("off")
    plt.savefig(output_file)
    plt.show()
    

if __name__ == '__main__':
    feature_column = int(sys.argv[1]) # zero based index
    count_column = int(sys.argv[2])
    output_file = sys.argv[3]
    input_files = sys.argv[4:]

    corpus = construct_corpus(input_files, feature_column, count_column)
    hdp = train_hdp(corpus, initial_k=5, iterations=2000)
    visualize_topics(get_topics(hdp))
    # hdp = LdaModel(corpus, id2word=dictionary, num_topics=5)
    # visualize_topics(hdp, 5)
    
