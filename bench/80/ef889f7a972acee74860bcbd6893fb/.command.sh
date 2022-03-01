#!/bin/bash -ue
cat C57BL_6NJ.H6.mm39.bed C57BL_6NJ.H7.mm39.bed | awk '!x[$0]++' > C57BL_6NJ.INS.validated.bed
