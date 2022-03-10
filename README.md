# Simon Lab Target Capture Pipeline

![Pipeline flowchart](/images/pipeline_flowchart.png)

This pipeline is for extracting target capture data from DNA assemblies and splitting paralogous loci. It incorporates [ALiBaSeq](https://github.com/AlexKnyshov/alibaseq). 

Pipeline steps:
...
Make BLAST databases of your reference genome(s).
Use `ref_genome_blast.sh` to find loci with paralogs that were missed by UPhO.
Use `paralog_remover.sh` to remove missed paralogs from the post-UPhO BLAST files.
Run ALiBaSeq using the post-UPhO BLAST files.
