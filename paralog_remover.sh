#Script that removes paralogs (sequences that did not overlap with same place in reference genome) from UPhO-split blast files
#Takes intersect_result_list.txt from ref_genome_blast.sh

#Usage message
if [[ $2 == "" ]]

then
  echo -e "To be run after ref_genome_blast.sh \nUsage: \
\n\$1 file containing list of paralogous sequences in loci \
\n\$2 folder containing UPhO-split blast files \nExample: bash paralog_remover.sh ./intersect_results/intersect_result_list.txt ./new_blastfiles"

else
  #Create output folder
  echo "Creating output folder \"blast_results_no_paralogs\""
  rm blast_results_no_paralogs/*
  mkdir blast_results_no_paralogs/

  #Removes each paralogous sequence from each blast file and saves file in output folder
  while read x
  do
    locus=$(echo $x | cut -f1 -d " ")
    taxon=$(echo $x | cut -f2 -d " ")
    echo "removing $locus from ${taxon}.blast"
    cp -n $2/${taxon}.blast ./blast_results_no_paralogs/
    sed -i '/${locus).fas/d' ./blast_results_no_paralogs/${taxon}.blast 
  done < $1

  #Copies remaining UPhO-split blast files to output folder
  echo "Copying remaining blast files to \"blast_results_no_paralogs\""
  cp -n $2/*.blast ./blast_results_no_paralogs/
fi
