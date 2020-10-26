#!/usr/bin/bash

#when Jason changes the ANN col from ; sep to \t sep, just take it out before you process it. You don't really need it anyway. 
n=$(head -n 1 Avian_CIFAR_5.snpEff.matrix.tsv | tr -dc \\t | tr \\t \\n | wc -l)
cut -f1-$n Avian_CIFAR_5.snpEff.matrix.tsv > Avian_CIFAR_5.snpEff.matrix.fixed.tab
