from collections import namedtuple
from gen_tree import main

if __name__ == "__main__":
    TestArgs = namedtuple('TestArgs', 'seq_file saito')
    args = TestArgs('./distance_utils/sequences.fasta', True)
    main(args)
    args = TestArgs('./distance_utils/sequences.fasta', False)
    main(args)
    