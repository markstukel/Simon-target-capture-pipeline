#Script to ensure post-UPhO locus alignments map to same region of genome
#Outputs list of all loci with at least one sequence that does not overlap
#Note: if a locus has a sequence that is in the reverse direction in the reference genome, bedtools will skip it with an error message.
#You may have to manually check skipped loci listed in stderr

#Usage statement
if [[ $3 == "" ]]

then
  echo -e "Usage: \n\$1 folder containing loci fasta files \n\$2 folder containing blast db(s) of reference genome(s) \n\$3 number of threads \
\nExample: bash ref_genome_blast.sh ./loci/ ./ref_genomes/ 16"

else
  #Create folder for blast results
  rm blast_results/*
  mkdir blast_results/

  #Make list of loci files for queries
  loci=$(ls $1/*.fas | rev | cut -f1 -d/ | rev | awk '$1=$1' ORS=' ')

  #Make list of reference genome blast dbs
  db=$(ls $2/*.fasta | rev | cut -f1 -d/ | rev | awk '$1=$1' ORS=' ')

  #Blast each locus query against each reference genome blast db
  for x in $loci; do 
    echo "Blasting $x";
    for y in $db; do 
      tblastx -query $1/$x -db $2/$y -outfmt 6 -max_target_seqs 1 -max_hsps 1 \
      -evalue 1e-5 -num_threads $3 -out blast_results/${x}_${y}.res;
      echo "  $y done";
    done
  done

  #Create folder for bedtools intersect results
  cd blast_results/
  
  rm intersect_results/*
  mkdir intersect_results/

  #Bedtools stuff
  echo "Finding remaining paralogs in loci"
  for result in *.res; do
    cat $result | awk '{print($2"\t"$9-1"\t"$10"\t"$1)}' | awk '{OFS="\t"} {if($2 > $3) {print $1, $3, $2, $4} else print $0 }' > intersect_results/${result}.bed; #Convert blast output to .bed file
    head intersect_results/${result}.bed -n 1 > intersect_results/${result}_B.bed; #Make comparison file for intersect from top line
    bedtools intersect -a intersect_results/${result}.bed -b intersect_results/${result}_B.bed \
    -loj > intersect_results/${result}_intersect.out; #Records overlap between each line of .bed file with first line and prints a null value if no overlap
  done
  cd intersect_results/
  
  #Makes list of sequences that mapped to different place in the reference than majority
  echo "Compiling paralogs in loci"
  for file in *.out
  do
    grep '\-1' $file > test1.txt
    grep -v '\-1' $file > test2.txt
    locus=$(echo $file | cut -f1 -d.)
    x=$((`wc -l < test1.txt`))
    y=$((`wc -l < test2.txt` / 2))
    if [[ $x -eq 0 ]]
    then
      :
    elif [[ $x -le $y ]]
    then
      while read line
        do
          taxon=$(echo $line | cut -f4 -d " " | cut -f1 -d\|)
          echo -e ${locus}"\t"${taxon} >> result_list.txt
        done < test1.txt
    else
      list2=$(grep -v '\-1' $file)
        while read line
        do
          taxon=$(echo $line | cut -f4 -d " " | cut -f1 -d\|)
          echo -e ${locus}"\t"${taxon} >> result_list.txt
        done < test2.txt
    fi
  done
  rm test1.txt
  rm test2.txt
  uniq result_list.txt > intersect_result_list.txt
  rm result_list.txt
fi
