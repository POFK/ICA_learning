#!/bin/bash
#git clone https://github.com/POFK/test.git
com=`date '+%Y/%m/%d_%H:%M:%S'`
git add .
git commit -m $com
#git commit -a 
git push origin master
