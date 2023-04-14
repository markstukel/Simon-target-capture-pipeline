#Script that uses UPHO to cluster and eliminate redundant sequences or close in-paralogs
#requires UPHO, seqkit
###Requires folder with final cdhit proccesed loci files, blast results folder, minimum number of taxa in a new cluster to include as new locus
###Outputs updated blast files in new_blastfiles folder for reprocessing and alignments post upho split by paralogs in post_upho_alignments

Upho_path=/home/FCAM/jvailionis/programs/UPhO/UPhO.py
module load seqkit


echo -e "Example: Upho.sh ~/Cicada_raw_reads/M_noflanks_out/aligned/finalcdhit /labs/Simon/AHE/Magicicada_AHE_processing/blast_results 10"
cd $1


echo -e "Running UPHO based on files in upho.files.txt "
list=$(cat upho.files.txt | awk '/>/ {printf("\n%s\n",$0);next; } { printf("%s"," "$0);}  END {printf("\n");}')
for x in $list
do
#mafft trees have _ in place of periods so must work around that 
sed -r -i 's/[[:digit:]]+_I/I/g' ../../$x'.fas.tree'
sed -i 's/_NODE/|NODE/g' ../../$x'.fas.tree'
sed -z -i "s/\n//g" ../../$x'.fas.tree'
python $Upho_path -in ../../$x'.fas.tree' -iP 
mv UPhO_orthogroups.csv $x.Upho.csv 
mv UPhO_nr_orthogroups.csv $x.Upho.nr.csv 
done


echo -e "Running UPHO parser"

#set up new folders
mkdir new_blastfiles
rm new_blastfiles/*
mkdir post_upho_alignments
rm post_upho_alignments/*


#Uphoparser script works on single locus given one cluster currently after fasta files have been renamed 
# it updates original blast files with new paralog assignment so alibaseq can retrieve original hits with flanks and parse out to paralogs 
##Upho parser 
for locus in $(ls L*.fas | sed 's/\..*//'); do
	if [ -f "$locus".Upho.csv ]; then
		if [ "$(wc -l <"$locus".Upho.csv)" -eq "$(wc -l <"$locus".Upho.nr.csv)" ]; then
			if [ "$(wc -l <"$locus".Upho.csv)" -eq 1 ]; then
				### UPHO GAVE ONE CLUSTER
				echo $locus upho, one cluster
				for blastfile in $2/*.fasta.blast; do
					outname="$(sed 's|.*/||' <<< $blastfile)"
					grep "$locus" "$blastfile" | grep -f <(awk '$1=$1' RS="," "$locus".Upho.nr.csv | tail -n +2 | cut -f 2 -d '|' | sed -r 's/(.*)_/\1./') >> new_blastfiles/$outname
					seqkit grep -f <(awk '$1=$1' RS="," "$locus".Upho.nr.csv | tail -n +2 | sed 's/_contigs_fasta/.contigs.fasta/' | sed -r 's/(.*)_/\1./') "$locus".fas > post_upho_alignments/"$locus".fas 
				done
			else
				#### UPHO GAVE MULTIPLE CLUSTERS
				echo $locus upho, multiple clusters
				while read -r line; do
					IFS=',' read -ra a <<< "$line"
					if [ "${#a[@]}" -gt $3 ]; then
						cluster_num=$(( $(echo "${a[0]}" | sed 's/#_//') +1 ))
						for seq in "${a[@]:1}"; do
							ass_name="$(cut -f 1 -d '|' <<< $seq | sed 's/_/./g')"
							node="$(cut -f 2 -d '|' <<< $seq | sed -r 's/(.*)_/\1./g')"
							#echo $locus $ass_name $node $cluster_num
							grep "$locus" $2/"$ass_name".blast | grep "$node" | sed 's/'"$locus"'/'"$locus"'_'"$cluster_num"'/' >> new_blastfiles/"$ass_name".blast
							seqkit grep -p "$ass_name"'|'"$node" "$locus".fas >> post_upho_alignments/"$locus"_"$cluster_num".fas
						done
					fi
				done < "$locus".Upho.nr.csv 
			fi
		else
			#### UPHO, 2 CLUSTERS IN REDUNDANT FILE BUT 1 NR CLUSTER IS WRONG 
			echo $locus upho, 2 r clusters but 1 in nr
			for blastfile in $2/*.fasta.blast; do
				outname="$(sed 's|.*/||' <<< $blastfile)"
				grep "$locus" "$blastfile" | grep -f <(head -n 1 "$locus".Upho.csv | awk '$1=$1' RS="," | tail -n +2 | cut -f 2 -d '|' | sed -r 's/(.*)_/\1./') >> new_blastfiles/$outname
				seqkit grep -f <(head -n 1 "$locus".Upho.csv | awk '$1=$1' RS="," | tail -n +2 | sed 's/_fasta/.fasta/' | sed -r 's/(.*)_/\1./') "$locus".fas > post_upho_alignments/"$locus".fas 
			done
		fi
	else
		#### NO UPHO ####
		echo $locus no upho
		targets="$(grep '>' "$locus".fas | tr -d '>')"
		for blastfile in $2/*.fasta.blast; do
			outname="$(sed 's|.*/||' <<< $blastfile)"
			grep "$locus" "$blastfile" | grep -f <(echo "$targets" | cut -f 2 -d "|") >> new_blastfiles/$outname
			cp "$locus.fas" post_upho_alignments/
		done
	fi
done
