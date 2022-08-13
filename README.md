SSHash-Lite
===========

A membership-only version of [SSHash](https://github.com/jermp/sshash). It works with files with duplicate k-mers too.

For instructions on how to compile the code,
build the dictionaries, and data format, please refer to [this](https://github.com/jermp/sshash) repository.

### Example usage

Build an index over the E. Coli matchtigs with

	./sshash-lite build ../data/ecoli_matchtigs_k31.fa.gz 31 15 -d tmp --check -o ecoli-matchtigs.sshash

and query with

	./sshash-lite query ecoli-matchtigs.sshash ../data/queries/SRR9873306_1.10K.fastq.gz -t 0.8

where a query read is considered as positive if at least 80% of the k-mers are positive.

### Enable verbose output

By default we print a summary report for the whole query file.
To print the result for each query, define the flag `SSHASH_QUERY_VERBOSE_OUTPUT ` as follows:

	cmake .. -D SSHASH_QUERY_VERBOSE_OUTPUT=On
	make -j
    
and re-run the query.

The query output is written on `std::cerr` by default, so it is possible to capture the output with

	./sshash-lite query ecoli-matchtigs.sshash ../data/queries/SRR9873306_1.10K.fastq.gz -t 0.8 2> query_output.txt
	
which writes the output to the file `query_output.txt`.