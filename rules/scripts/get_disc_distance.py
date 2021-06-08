import sys
import numpy as np

insert_sizes = [float(i.rstrip()) for i in open(sys.argv[1], 'r').readlines()]
mean = np.mean(insert_sizes)
std = np.std(insert_sizes)

discordant_distance = int(mean + 4*std)
print(discordant_distance, end='')
