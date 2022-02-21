# echo "new" > snapshots.txt
echo "genome ${igv_workdir}/Mus_musculus.GRCm39.dna.toplevel.fa" >> snapshots.txt
echo "snapshotDirectory ./${dir}" >> snapshots.txt

echo "load ${igv_workdir}/${strain}.session.xml" >> snapshots.txt
echo "load ${bed_file}" >> snapshots.txt

bedToIgv -slop 50 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt
bedToIgv -slop 500 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt

echo "exit" >> snapshots.txt

sed -i -e "s/.png/_-_${simple_name}.png/g" snapshots.txt
sed -i -e "s/snapshot /snapshot ${type}/g" snapshots.txt

echo "IGV.Bounds=0,0,2304,1296" > prefs.properties
# echo "SAM.SHOW_ALL_BASES=false" >> prefs.properties
# echo "SAM.SHOW_SOFT_CLIPPED=true" >> prefs.properties
# echo "SAM.SHOW_JUNCTION_TRACK=false" >> prefs.properties
# echo "SAM.SHOW_JUNCTION_FLANKINGREGIONS=false" >> prefs.properties
echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties


xvfb-run --auto-servernum -s "-screen 0 2304x1296x24" java -Xmx2000m --module-path=/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties