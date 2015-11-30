#!/usr/bin/env python
# coding=utf-8
from subprocess import call
import os
list=os.listdir('./')
list.sort()
call('cp ./result.tex.model result.tex',shell=True)
f1 = open('./result.tex', 'a')
name=[]
f1.writelines('\\begin{center}')
for i in list:
    if i[-4:]=='.eps':
        f1.writelines('\includegraphics[scale=0.39]{./'+i+'}'+'\n')
        name.append(i)
f1.writelines('\end{center}')
f1.writelines('\end{document}'+'\n')
f1.close()
call('pdflatex result.tex',shell=True)
call('rm *aux *log *.out *to.pdf',shell=True)
call('evince result.pdf',shell=True)
