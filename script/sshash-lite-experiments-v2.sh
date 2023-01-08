# !/usr/bin/sh

#### Human dataset

# regular index
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_greedytigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-greedytigs.sshash --verbose > results-28-10-22/human-greedytigs.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_unitigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-unitigs.sshash --verbose > results-28-10-22/human-unitigs.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_simplitigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-simplitigs.sshash --verbose > results-28-10-22/human-simplitigs.indexing_log

./sshash-lite query human-greedytigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-greedytigs.query_log
./sshash-lite query human-greedytigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-greedytigs.query_log
./sshash-lite query human-unitigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-unitigs.query_log
./sshash-lite query human-unitigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-unitigs.query_log
./sshash-lite query human-simplitigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-simplitigs.query_log
./sshash-lite query human-simplitigs.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-simplitigs.query_log

# canonical index
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_greedytigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-greedytigs-canon.sshash --verbose > results-28-10-22/human-greedytigs-canon.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_unitigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-unitigs-canon.sshash --verbose > results-28-10-22/human-unitigs-canon.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/human_pangenome_simplitigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-simplitigs-canon.sshash --verbose > results-28-10-22/human-simplitigs-canon.indexing_log

./sshash-lite query human-greedytigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-greedytigs-canon.query_log
./sshash-lite query human-greedytigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-greedytigs-canon.query_log
./sshash-lite query human-unitigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-unitigs-canon.query_log
./sshash-lite query human-unitigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-unitigs-canon.query_log
./sshash-lite query human-simplitigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 > results-28-10-22/human-simplitigs-canon.query_log
./sshash-lite query human-simplitigs-canon.sshash /data2/DNA/human_salmonella300K/human_pangenome_query.fa.gz -t 0.8 >> results-28-10-22/human-simplitigs-canon.query_log



#### Salmonella-300K dataset

# regular index
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_greedytigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-greedytigs.sshash --verbose > results-28-10-22/salmonella300k-greedytigs.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_unitigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-unitigs.sshash --verbose > results-28-10-22/salmonella300k-unitigs.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_simplitigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-simplitigs.sshash --verbose > results-28-10-22/salmonella300k-simplitigs.indexing_log

./sshash-lite query salmonella300k-greedytigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-greedytigs.query_log
./sshash-lite query salmonella300k-greedytigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-greedytigs.query_log
./sshash-lite query salmonella300k-unitigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-unitigs.query_log
./sshash-lite query salmonella300k-unitigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-unitigs.query_log
./sshash-lite query salmonella300k-simplitigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-simplitigs.query_log
./sshash-lite query salmonella300k-simplitigs.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-simplitigs.query_log

# canonical index
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_greedytigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-greedytigs-canon.sshash --verbose > results-28-10-22/salmonella300k-greedytigs-canon.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_unitigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-unitigs-canon.sshash --verbose > results-28-10-22/salmonella300k-unitigs-canon.indexing_log
./sshash-lite build /data2/DNA/human_salmonella300K/salmonella_simplitigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-simplitigs-canon.sshash --verbose > results-28-10-22/salmonella300k-simplitigs-canon.indexing_log

./sshash-lite query salmonella300k-greedytigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-greedytigs-canon.query_log
./sshash-lite query salmonella300k-greedytigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-greedytigs-canon.query_log
./sshash-lite query salmonella300k-unitigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-unitigs-canon.query_log
./sshash-lite query salmonella300k-unitigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-unitigs-canon.query_log
./sshash-lite query salmonella300k-simplitigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 > results-28-10-22/salmonella300k-simplitigs-canon.query_log
./sshash-lite query salmonella300k-simplitigs-canon.sshash /data2/DNA/human_salmonella300K/salmonella_query.fa.gz -t 0.8 >> results-28-10-22/salmonella300k-simplitigs-canon.query_log