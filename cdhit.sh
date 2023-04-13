#Script that uses cdhit to cluster and eliminate redundant sequences or close in-paralogs
#requires cdhit/4.6.8
#Takes aligned folder of no flanks output of many .fas files 
#Outputs locus file with removed contigs as Lx.merged.fas within cdhit# folder then overwrites new files and combines with old files into finalcdhit folder
#Internal adjustable defaults cluster threshold at 97% 
#clt=.97
echo -e "Example: cdhit.sh ./M_noflanks_out/aligned"


#sort out dupes
# works for non merging single duplicate....works for non merging double duplicates...works for merging one duplicate
# merge files will be empty if nothing accomplished
# make multiline fasta into single line fasta for each file
echo -e "Converting multiline fasta files into single line temporary fasta files" 
cd $1
for x in *.fas
do
#replace end of line to string concurrent sequence lines together
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $x > "$x".oneline
done


echo -e "Running CDHIT on all contigs with multiple representatives from one taxon per alignment " 
#stripped clt
clt2=$(echo $clt | cut  -d. -f 2)
module load cdhit/4.6.8
#loop across all loci
for x in *.fas
do
#add taxa with duplicates to temporary taxa file
grep ">" $x | cut -f 1 -d "|" | sort | uniq -d > "x".taxa
#read temporary taxa file and perform cd-hit on all sequences from the same taxon
IFS=$'\n' read  -d'' -r -a taxa  < "x".taxa
for i in "${taxa[@]}"
do
san=$(echo $i| cut -f 2 -d ">")
grep -A 1 -h $san "$x".oneline > "$x"."$san".temp
cd-hit-est -i "$x"."$san".temp -o "$x"."$san".cdhit -c $clt -G 0 -A 15 -aS .01 -aL 0.01 -d 80
done
done
rm x.taxa

echo -e "Moving cdhit clstr files to subfolder " 
mkdir cdhit$clt2
mv *.cdhit.clstr ./cdhit$clt2/

echo -e "Splicing remaining contigs into merged files" 
for x in *.fas
do
w=$( echo $x | cut -f 1 -d .)
#Find secondary contigs to get remove per taxon per locus (those consisting of the second and shorter sequence that is above the threshold of matching)
rm=$(grep -E -v  '^>|^0' ./cdhit$clt2/"$x"*.cdhit.clstr | cut -f 2 -d " " | sed 's/\.\.\.//g' - | cut -f 2 -d "|")
rm2=$(echo $rm | sed  's/ /|/g')
#Delete secondary contigs and print remaining sequences to new merged file
awk '/('$rm2')/{getline;next}{print}' $x.oneline > ./cdhit$clt2/$w'.merged.fas'
done

echo -e "Cleaning temporary files"
rm *oneline
rm *.temp
rm *.cdhit
#remove empty merged files
find ./cdhit$clt2/ -size  0 -print -delete

 
#count how many merged 
f=$(grep ^1 ./cdhit$clt2/*cdhit.clstr | wc -l)
echo -e $f "files with some secondary contigs more than" $clt2 "% identical removed"

echo -e "condensing merged files and old files into finalcdhit folder"
#Overwrite old files with new merged files into processed folder and then into new folder for next step 
cd  ./cdhit$clt2/
mkdir ../processedcdhit$clt2/
for x in *.merged.fas
do 
y=$(echo $x | cut -f 1 -d .)
#echo $y
cp $x ../processedcdhit$clt2/$y.fas
done
cd ..
mkdir finalcdhit
cp *.fas ./finalcdhit/
cp ./processedcdhit$clt2/*.fas ./finalcdhit/
cd ./finalcdhit

echo -e "checking list of loci with mutliple contigs per taxon to process with UpHO"
# look for files retaining duplicates still 
for x in *.fas; do echo $x; grep ">" $x | cut -f 1 -d "|" | sort | uniq -d | grep -B 1 ">" | grep -r ^L -; done > upho.files.txt
#Remove .fas extension from upho.files.txt
sed -i 's/.fas//g' upho.files.txt
z=$(wc -l upho.files.txt)
echo -e $z "remaining loci with mutliple contigs per taxon to process with UpHO"
