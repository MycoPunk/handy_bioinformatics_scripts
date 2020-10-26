#!/usr/bin/bash

n=$(head -n 1 Avian_CIFAR_5.snpEff.matrix.tsv | tr -dc \\t | tr \\t \\n | wc -l)
cut -f1-$n Avian_CIFAR_5.snpEff.matrix.tsv > Avian_CIFAR_5.snpEff.matrix.fixed.tab
