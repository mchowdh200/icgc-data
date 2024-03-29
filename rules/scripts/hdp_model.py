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

def train_hdp(
        corpus,
        initial_k=10,
        term_weight=tp.TermWeight.PMI,
        gamma=1, alpha=0.1,
        iterations=2000):
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
        plt.subplot(nrows, ncols, i+1)
        freqs = {word : freq for word, freq in topic}
        wc = WordCloud(background_color="black", max_words=100)
        wc.generate_from_frequencies(freqs)
        plt.imshow(wc, interpolation="bilinear")
        plt.axis("off")
        # plt.barh(range(len(freqs)), freqs.values())
        # plt.yticks(range(len(freqs)), list(freqs.keys()))
    plt.savefig(output_file)
    plt.show()

def eval_coherence(model, metric='c_v', top_n=25):
    return tp.coherence.Coherence(
        model, coherence=metric, top_n=top_n).get_score()



if __name__ == '__main__':
    feature_column = int(sys.argv[1]) # zero based index
    count_column = int(sys.argv[2])
    output_file = sys.argv[3]
    input_files = sys.argv[4:]

    corpus = construct_corpus(input_files, feature_column, count_column)

    ## grid search
    best_model = None
    best_coherence = 0
    for term_weight in [tp.TermWeight.ONE,
                        tp.TermWeight.PMI,
                        tp.TermWeight.IDF]:
        for alpha in [0.01, 0.1, 0.5]:
            for gamma in [0.1, 0.5, 1.0]:
                print('TRAINING MODEL WITH PARAMS')
                print(f'{term_weight = }\t{alpha = }\t{gamma = }')
                model = train_hdp(
                    corpus, initial_k=5, term_weight=term_weight,
                    gamma=gamma, alpha=alpha, iterations=2000)

                coherence = eval_coherence(model)
                print(f'{coherence = }')
                if best_model:
                    if coherence > best_coherence:
                        best_model = model
                        best_coherence = coherence
                        best_weight = term_weight
                        best_alpha = alpha
                        best_gamma = gamma
                else:
                    best_model = model
                    best_coherence = coherence
                    continue
                print('BEST MODEL PARAMS')
                print(f'{best_weight = }\t{best_alpha = }\t{best_gamma = }')
                print(f'{best_coherence = }')
                        
    print('SAVING MODEL')
    best_model.save('hdp_model.bin')
    visualize_topics(get_topics(best_model))
