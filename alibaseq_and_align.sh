#Script that uses alibaseq tool to parse blast files and pull out matching contigs twice with flanks included the second time. Contigs without flanks are aligned for use in orthology refinement.
#Requires alibaseq path and MAFFT
#Takes blast file folder, target assembly folder and name
#Internal adjustable defaults are blast evalue of 1e15, flanks of size 300 bp and using 32 threads for MAFFT alignment

#Usage message
if [[ $3 == "" ]]

then
  echo -e "To be run after folder_blast.sh \nUsage: \
\n\$1 folder containing blast results \
\n\$2 folder containing assemblies \
\n\$3 name of dataset \nExample: bash alibaseq_and_align.sh /labs/Simon/AHE/Magicicada_AHE_processing/blast_results /labs/Simon/AHE/Magicicada_AHE_processing/Magicicada_assemblies Mag"



else
  echo "Running no flanks alibaseq"
  python /home/FCAM/egordon/alibaseq/alibaseq.py -b $1 -f M -x s -t $2 -e 1e-15 --ac aa-tdna -c 0
  echo "Creating no flanks output folder \"$3_no_flanks_out\" and moving output and log files"
  mv alibaseq_out $3_noflanks_out
  mkdir $3_noflanks_out/logs
  mv *default* $3_noflanks_out/logs
  echo "Running 300 bp flanks alibaseq"
  python /home/FCAM/egordon/alibaseq/alibaseq.py -b $1 -f M -x s --fl 300 -t $2 -e 1e-15 --ac aa-tdna -c 0
  echo "Creating 300 bp flanks output folder \"$3_300bpflanks_out\" and moving output and log files"
  mv alibaseq_out $3_300bpflanks_out
  mkdir $3_300bpflanks_out/logs
  mv *default* $3_300bpflanks_out/logs
  echo "Aligning no flanks alibaseq output in subdirectory aligned folder with 32 threads"
  cd $3_noflanks_out
  mkdir aligned
  for x in *.fas
        do
	einsi --treeout --thread 32  --averagelinkage $x > ./aligned/$x
        done
fi
