#!/bin/bash -ue
cat C57BL_6NJ.H1.mm39.bed C57BL_6NJ.H2.mm39.bed | awk '!x[$0]++' > C57BL_6NJ.DEL.validated.bed
