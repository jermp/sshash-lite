# !/usr/bin/sh

#### Human dataset

# regular index
./sshash-lite build /data2/DNA/human/greedytigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-greedytigs.sshash --verbose > human-greedytigs.indexing_log
./sshash-lite build /data2/DNA/human/unitigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-unitigs.sshash --verbose > human-unitigs.indexing_log
./sshash-lite build /data2/DNA/human/usttigs.fa.gz 31 20 -d /data2/DNA/tmp_dir -o human-usttigs.sshash --verbose > human-usttigs.indexing_log

./sshash-lite query human-greedytigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-greedytigs.query_log
./sshash-lite query human-greedytigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-greedytigs.query_log
./sshash-lite query human-unitigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-unitigs.query_log
./sshash-lite query human-unitigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-unitigs.query_log
./sshash-lite query human-usttigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-usttigs.query_log
./sshash-lite query human-usttigs.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-usttigs.query_log

# canonical index
./sshash-lite build /data2/DNA/human/greedytigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-greedytigs-canon.sshash --verbose > human-greedytigs-canon.indexing_log
./sshash-lite build /data2/DNA/human/unitigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-unitigs-canon.sshash --verbose > human-unitigs-canon.indexing_log
./sshash-lite build /data2/DNA/human/usttigs.fa.gz 31 19 -d /data2/DNA/tmp_dir --canonical-parsing -o human-usttigs-canon.sshash --verbose > human-usttigs-canon.indexing_log

./sshash-lite query human-greedytigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-greedytigs-canon.query_log
./sshash-lite query human-greedytigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-greedytigs-canon.query_log
./sshash-lite query human-unitigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-unitigs-canon.query_log
./sshash-lite query human-unitigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-unitigs-canon.query_log
./sshash-lite query human-usttigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 > human-usttigs-canon.query_log
./sshash-lite query human-usttigs-canon.sshash /data2/DNA/human/query.fa.gz -t 0.8 >> human-usttigs-canon.query_log



#### Salmonella-300K dataset

# regular index
./sshash-lite build /data2/DNA/salmonella300k/greedytigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-greedytigs.sshash --verbose > salmonella300k-greedytigs.indexing_log
./sshash-lite build /data2/DNA/salmonella300k/unitigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-unitigs.sshash --verbose > salmonella300k-unitigs.indexing_log
./sshash-lite build /data2/DNA/salmonella300k/usttigs.fa.gz 31 17 -d /data2/DNA/tmp_dir -o salmonella300k-usttigs.sshash --verbose > salmonella300k-usttigs.indexing_log

./sshash-lite query salmonella300k-greedytigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-greedytigs.query_log
./sshash-lite query salmonella300k-greedytigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-greedytigs.query_log
./sshash-lite query salmonella300k-unitigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-unitigs.query_log
./sshash-lite query salmonella300k-unitigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-unitigs.query_log
./sshash-lite query salmonella300k-usttigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-usttigs.query_log
./sshash-lite query salmonella300k-usttigs.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-usttigs.query_log

# canonical index
./sshash-lite build /data2/DNA/salmonella300k/greedytigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-greedytigs-canon.sshash --verbose > salmonella300k-greedytigs-canon.indexing_log
./sshash-lite build /data2/DNA/salmonella300k/unitigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-unitigs-canon.sshash --verbose > salmonella300k-unitigs-canon.indexing_log
./sshash-lite build /data2/DNA/salmonella300k/usttigs.fa.gz 31 16 -d /data2/DNA/tmp_dir --canonical-parsing -o salmonella300k-usttigs-canon.sshash --verbose > salmonella300k-usttigs-canon.indexing_log

./sshash-lite query salmonella300k-greedytigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-greedytigs-canon.query_log
./sshash-lite query salmonella300k-greedytigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-greedytigs-canon.query_log
./sshash-lite query salmonella300k-unitigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-unitigs-canon.query_log
./sshash-lite query salmonella300k-unitigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-unitigs-canon.query_log
./sshash-lite query salmonella300k-usttigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 > salmonella300k-usttigs-canon.query_log
./sshash-lite query salmonella300k-usttigs-canon.sshash /data2/DNA/salmonella300k/query.fa.gz -t 0.8 >> salmonella300k-usttigs-canon.query_log