#Script that uses cdhit to cluster and eliminate redundant sequences or close in-paralogs
#requires seqkit and cdhit/4.6.8
#Takes aligned folder of no flanks output of many .fas files 
#Outputs locus file with removed contigs as Lx.merged.fas
#Internal adjustable defaults cluster threshold at 97% 
#clt=.97
echo -e "Example: cdhit.sh ./M_noflanks_out/aligned"


#sort out dupes
# works for non merging single duplicate....works for non merging double duplicates...works for merging one duplicate
# merge files will be empty if nothing accomplished
# make multiline fasta into single line fasta for each file
echo -e "Converting multiline fasta files into single line temporary fasta files" 
for x in *.fas
do
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $x > "$x".oneline
done


echo -e "Running CDHIT on all contigs with multiple representatives from one taxon per alignment " 
#stripped clt
clt2=$(echo $clt | cut  -d. -f 2)
for x in *.fas
do
#add taxa with duplicates to taxa array
grep ">" $x | cut -f 1 -d "|" | sort | uniq -d > "x".taxa
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
#contigs to get rid of 
rm=$(grep ^1 ./cdhit$clt2/"$x"*.cdhit.clstr | cut -f 2 -d " " | sed 's/\.\.\.//g' - | cut -f 2 -d "|")
rm2=$(echo $rm | sed  's/ /|/g')
awk '/('$rm2')/{getline;next}{print}' $x.oneline > ./cdhit$clt2/$w'.merged.fas'
done

echo -e "Cleaning temporary files"
rm *oneline
rm *.temp
rm *.cdhit
#count how many merged ......84 with 100%, 331 with 99%, 429 with 98%, 512 with 97% , 603 with 95%, 809 with 90%...817 with 80  5083 total
grep ^1 *cdhit.clstr | wc
#remove empty merged files
find ./cdhit$clt2/ -size  0 -print -delete

 
#count how many merged 
f=$(grep ^1 ./cdhit$clt2/*cdhit.clstr | wc -l)
echo -e $f "files with some contigs removed"

# look for files retaining duplicates still 
for x in *.fas; do echo $x; grep ">" $x | cut -f 1 -d "|" | sort | uniq -d ; done > dupcount.txt
grep -B 1 ">" dupcount.txt | grep -r ^L - > upho.files.txt

