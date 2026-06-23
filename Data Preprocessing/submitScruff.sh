batch=$1
file=$2

# Because I create the full path to the input file in ./runScruff.sh dynamically using $batch, 
# I can put the relative path as the initial argument.
qsub -P pcga -N Scruff_Batch$batch -o ./scruff_log/Batch${batch}.log -j y -v BATCH=$batch -v FILE=$file ./runScruff.sh
