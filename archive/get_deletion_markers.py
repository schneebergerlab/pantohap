from util import updatesnpgenotypes
import sys
"""
Required arguments:
pwd: Working directory
c: chromosome ID
nproc: Number of processes
"""
updatesnpgenotypes(*sys.argv[1:])
