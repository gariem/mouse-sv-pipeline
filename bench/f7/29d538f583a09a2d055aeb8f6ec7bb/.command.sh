#!/bin/bash -ue
cat C3H_HeJ.H1.mm39.bed C3H_HeJ.H2.mm39.bed | awk '!x[$0]++' > C3H_HeJ.DEL.validated.bed
