#!/usr/bin/bash
pigz -cd -p 2 $1 | sed -n 2~4p |cut -c3-27 | paste - <(pigz -cd -p 2 $2 | sed 3~4d | sed 'N;N;s/\n/\t/g'  ) |  perl $3/trim.pl 
