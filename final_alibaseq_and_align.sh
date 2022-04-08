#Script that uses alibaseq tool to parse curated blast files and pull out matching contigs with flanks. Contigs are aligned with MAFFT.
#Requires alibaseq path and MAFFT
#Takes curated blast file folder, target assembly folder and name
#Internal adjustable defaults are blast evalue of 1e15, flanks of size 300 bp and using 32 threads for MAFFT alignment


##Defaults:
#Set evlaue
e='1e-15'
# set flank size
f='300'
# set thread number
t='32'

#alibaseq path
ali_path='/home/FCAM/egordon/alibaseq/alibaseq.py'


#Usage message
if [[ $3 == "" ]]

then
  echo -e "To be run after paralog_remover.sh \nUsage: \
\n\$1 folder containing blast results \
\n\$2 folder containing assemblies \
\n\$3 name of dataset \nExample: bash alibaseq_and_align.sh ~/Cicada_raw_reads/M_noflanks_out/aligned/finalcdhit/blast_results_no_paralogs/ /labs/Simon/AHE/Magicicada_AHE_processing/Magicicada_assemblies Mag"



else
  echo "Running no flanks alibaseq"
  python $ali_path -b $1 -f M -x s -t $2 -e $e --ac aa-tdna -c 0

  mkdir $3_noflanks_out/logs
  mv *default* $3_noflanks_out/logs
  echo "Running" $f "bp flanks alibaseq"
  python $ali_path -b $1 -f M -x s --fl $f -t $2 -e $e --ac aa-tdna -c 0
  echo "Creating "$f "bp flanks output folder \"$3_$fbpflanks_out\" and moving output and log files"
  mv alibaseq_out $3_$fbpflanks_out
  mkdir $3_$fbpflanks_out/logs
  mv *default* $3_$fbpflanks_out/logs
  echo "Aligning no flanks alibaseq output in subdirectory aligned folder with" $t "threads"
  cd $3_noflanks_out
  mkdir aligned
  for x in *.fas
        do
	einsi --treeout --thread $t  --averagelinkage $x > ./aligned/$x
        done
fi
