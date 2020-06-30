gzip -d merge.fna.gz
iss generate --cpus 4 -g merge.fna -b abundance_file.txt --n_reads 5M -m Hiseq --output test_data
