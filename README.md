# Simon Lab Target Capture Pipeline

![Pipeline flowchart](/images/pipeline_flowchart.png)

This pipeline is for extracting target capture data from DNA assemblies and splitting paralogous loci. It incorporates [ALiBaSeq](https://github.com/AlexKnyshov/alibaseq). 

Pipeline steps:

...

Use `folder_blast.sh` to Make BLAST databases of your reference genome(s) and search them. (?)

Use `alibaseq_and_align.sh` to parse blast files and pull out matching contigs twice (with flanks included the second time). Contigs without flanks are aligned for upstream use in orthology refinement.

Use `CDHIT.sh` to cluster and eliminate redundant sequences or close in-paralogs.

Use `UPhO.sh` to cluster into putative ortholog groups and eliminate paralogs updating blast files.

Use `ref_genome_blast.sh` to find loci with paralogs that were missed by UPhO.

Use `paralog_remover.sh` to remove missed paralogs from the post-UPhO BLAST files.

Use `final_alibaseq_and_align.sh` to run ALiBaSeq using the post-UPhO BLAST files.
