#!/bin/bash -ue
cat C3H_HeJ.H6.mm39.bed C3H_HeJ.H7.mm39.bed | awk '!x[$0]++' > C3H_HeJ.INS.validated.bed
