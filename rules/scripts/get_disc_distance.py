import sys
import numpy as np

insert_sizes = [float(i.rstrip()) for i in open(sys.argv[1], 'r').readlines()]
median = np.median(insert_sizes)
median_absolute_deviation = np.median(np.abs(insert_sizes - median))
discordant_distance = int(median + 4*median_absolute_deviation)
print(discordant_distance, end='')
