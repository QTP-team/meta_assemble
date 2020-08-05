gzip -d merge_sp3.fna.gz
iss generate --compress --cpus 8 -g merge_sp3.fna -b abundance_file.txt --n_reads 5M -m Hiseq --output test_data
rm *fastq
