conda activate picrust2

picrust2_pipeline.py -s otu.fasta -i otu_table.txt -o picrust2 -p 48

zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
python3 ${db} /home/xiaoxiutan/picrust/EasyMicrobiome/script/summarizeAbundance.py -i KEGG.KO.txt -m ${db} /home/xiaoxiutan/picrust/EasyMicrobiome/kegg/KO1-4.txt -c 2,3,4 -s ',+,+,' -n raw -o KEGG

python3 ${db} script/summarizeAbundance.py -i KEGG.KO.txt -m ${db} kegg/KO1-4.txt -c 2,3,4 -s ',+,+,' -n raw -o KEGG

