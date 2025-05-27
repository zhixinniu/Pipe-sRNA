# miRNA QC
# Software: miRTrace
# Installation: conda install mirtrace
# https://github.com/friedlanderlab/mirtrace
java -Xms4G -Xmx4G -jar /PATH/TO/PKG/mirtrace.jar qc -t 30 -s hsa *.gz  # run qc mode
java -Xms4G -Xmx4G -jar /PATH/TO/PKG/mirtrace.jar trace -t 30 -s hsa *.gz  # run trace mode
# or
mirtrace qc -t 30 -s hsa *.gz  # run qc mode
mirtrace trace qc -t 30 -s hsa *.gz  # run trace mode

# Manual: https://github.com/friedlanderlab/mirtrace/blob/master/release-bundle-includes/doc/manual/mirtrace_manual.pdf