#!/usr/bin/bash

paste <(gzip -cd $(echo $1 | tr , ' ') | sed -n '1~4p;2~4p' ) <(gzip -cd $(echo $2 | tr , ' ') | sed -n '2~4p;4~4p')  | sed 'N;s/\n/\t/' | perl $3/prepare_bam.pl 


