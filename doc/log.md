## Debug log

Jul. 2026

- The neighbor finder is now checked to be correct. OpenMP is easier in c++ because the variables can be defined within a loop, the advantage of `threads` lib is that no extra-lib is needed. However, it is still a good way to practice because if the code is runable using the `threads`, they can easily run using OpenMP. While if I do not thinking about the `threads` while coding, the code will probably not run using OpenMP. So I dicide to keep on the threads.

- About rdf, check the histogram constructing part.

- Memory error in NF is because the set pbc part did not set correct bins as boundary bin, so the probing bin is out of the box.

## Parallel note

