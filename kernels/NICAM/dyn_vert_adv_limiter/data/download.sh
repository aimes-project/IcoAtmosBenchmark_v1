#! /bin/bash

# please get from http://scale.aics.riken.jp/yashiro/IAB_database/snapshot.dyn_vert_adv_limiter.pe000000

file=snapshot.dyn_vert_adv_limiter.pe000000

rm -f ./$file
wget http://scale.aics.riken.jp/yashiro/IAB_database/$file

echo "Checking md5sum:"
md5sum -c ${file}.md5
