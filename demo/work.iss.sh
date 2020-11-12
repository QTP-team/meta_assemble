gzip -d merge_sp3.fna.gz
iss generate --compress --cpus 8 -g merge_sp3.fna -b abundance1_file.txt --n_reads 5M -m Hiseq --output test1_data
iss generate --compress --cpus 8 -g merge_sp3.fna -b abundance2_file.txt --n_reads 5M -m Hiseq --output test2_data
iss generate --compress --cpus 8 -g merge_sp3.fna -b abundance3_file.txt --n_reads 5M -m Hiseq --output test3_data
rm *fastq
