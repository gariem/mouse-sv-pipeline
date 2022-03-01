#!/bin/bash -ue
cat A_J.H2.mm39.bed A_J.H1.mm39.bed | awk '!x[$0]++' > A_J.DEL.validated.bed
