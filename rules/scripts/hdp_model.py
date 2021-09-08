import sys
import numpy as np
from pprint import pprint
from gensim.models import HdpModel, LdaModel
from gensim.corpora import Dictionary
import matplotlib.pyplot as plt
from wordcloud import WordCloud

def construct_corpus(input_files: list[str], feature_column: int,
                     count_column: int) -> tuple[list[list[str]], Dictionary]:
    """
    Contstruct a list of samples where each sample will contain
    a list of features present in that sample (stored as strings).
    Returns 
    """
    samples = []
    for bed in input_files:
        S = []
        with open(bed) as f:
            for line in f:
                if int((x:=line.rstrip().split())[count_column]) > 0:
                    S.extend([x[feature_column]] * int(x[count_column]))
            samples.append(S)
                
    # samples = [[x[feature_column]
    #             for line in open(bed).readlines()
    #             if int((x:=line.rstrip().split())[count_column]) > 0]
    #            for bed in input_files]

    dictionary = Dictionary(samples)
    dictionary.filter_extremes(no_below=2, no_above=0.5)
    bag_of_words = [dictionary.doc2bow(s) for s in samples]
    return bag_of_words, dictionary

def visualize_topics(topic_model: HdpModel, num_topics: int) -> None:
    """
    Generate word cloud visualizations for each topic in
    the provided topic model (saved as png).
    """
    topics = topic_model.show_topics(
        num_topics=num_topics,
        num_words=100,
        formatted=False)

    # list of topic word weights for each topic
    weights : list[dict[string, float]] = []
    for topic in topics:
        words = topic[1] # list of (words, weights) in topic
        weights.append({word:weight for word, weight in words})
    
    nrows = 2
    ncols = np.ceil(num_topics/2).astype(int)
    for i, topic in enumerate(weights):
        wc = WordCloud(background_color="black", max_words=1000)
        wc.generate_from_frequencies(topic)
        plt.subplot(nrows, ncols, i+1)
        plt.imshow(wc, interpolation="bilinear")
        plt.axis("off")
    plt.savefig(output_file)
    plt.show()
    

if __name__ == '__main__':
    feature_column = int(sys.argv[1]) # zero based index
    count_column = int(sys.argv[2])
    output_file = sys.argv[3]
    input_files = sys.argv[4:]

    corpus, dictionary = construct_corpus(input_files,
                                          feature_column,
                                          count_column)
    # print(corpus)
    # hdp = HdpModel(corpus, dictionary)
    # visualize_topics(hdp, hdp.m_T)
    hdp = LdaModel(corpus, id2word=dictionary, num_topics=10)
    visualize_topics(hdp, 10)
    
