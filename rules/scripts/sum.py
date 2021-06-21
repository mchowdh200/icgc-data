import sys

print(sum(int(i) for line in sys.stdin
          for i in line.rstrip().split()))
