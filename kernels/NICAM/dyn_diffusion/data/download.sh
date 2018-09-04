#! /bin/bash

# please get from http://r-ccs-climate.riken.jp/members/yashiro/IAB_database/snapshot.dyn_diffusion.pe000000

file=snapshot.dyn_diffusion.pe000000

rm -f ./$file
wget http://r-ccs-climate.riken.jp/members/yashiro/IAB_database/$file || exit 1

echo "Checking md5sum:"
md5sum -c ${file}.md5 || exit 1
