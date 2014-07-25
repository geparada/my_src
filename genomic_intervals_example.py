from random import randint, seed

# if you can install bx python then uncomment the line below
#
#from bx.intervals.operations.quicksect import IntervalNode

# otherwise just download the quickset module as shown above 
# and place it in next to your program
#
from quicksect import IntervalNode

# the span of the generated intervals
SPAN = 10

# the size of the genome
SIZE = 5*10**4

# the number of intervals
N = 10**4

def generate(x):
    "Generates random interval over a size and span"
    lo = randint(10000, SIZE)
    hi = lo + randint(1, SPAN)
    return (lo, hi)

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]

# use this to force both examples to generate the same data
seed(10)

# generate 10 thousand random intervals
data = map(generate, xrange(N))

# generate the intervals to query over
query = map(generate, xrange(10))

# start the root at the first element
start, end = data[0]
tree = IntervalNode( start, end )

# build an interval tree from the rest of the data
for start, end in data[1:]:
    tree = tree.insert( start, end )

for start, end in query:
    overlap = find(start, end , tree)
    print '(%s, %s) -> %s' % (start, end, overlap)