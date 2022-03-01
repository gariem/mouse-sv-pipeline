#!/bin/bash -ue
cat A_J.H6.mm39.bed A_J.H7.mm39.bed | awk '!x[$0]++' > A_J.INS.validated.bed
