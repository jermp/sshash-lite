import sys, gzip as gz
import io

input_multifasta_filename = sys.argv[1]

with io.TextIOWrapper(gz.open(input_multifasta_filename)) as f:
    res = ""
    for line in f:
        l = line.strip()
        if ">" in l:
            if res:
                print(res)
            print(line, end='')
            res = ""
        else:
            res += l
    print(res)