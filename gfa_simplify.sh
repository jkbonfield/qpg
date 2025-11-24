#!/bin/sh

file=$1
M=${2:-5}
P=`awk '/^P/ {print $2;exit}' $file`
echo Running "vg simplify -P $P -m $M $file | odgi unchop | odgi view" 1>&2
vg simplify -m $M -P $P $file > $file.simplify.gfa
odgi unchop -i $file.simplify.gfa -o $file.simplify.og
odgi view -g -i $file.simplify.og > $file.simplify.gfa

# pathfinder can't handle L records before S records.
egrep '^[HS]' $file.simplify.gfa
egrep '^L'    $file.simplify.gfa

odgi stats -S -L -i $file.simplify.og 1>&2



# See also:
#smoothxg -M -i $file -o $file.smooth -r40 --consensus-spec=$file,100
#gfaffix


