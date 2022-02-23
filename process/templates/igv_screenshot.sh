echo "snapshotDirectory ./${dir}" >> snapshots.txt

echo "load ${igv_workdir}/${strain}/${strain}.session.xml" >> snapshots.txt
echo "load ${bed_file}" >> snapshots.txt

bedToIgv -slop 50 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt
bedToIgv -slop 500 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.txt

echo "exit" >> snapshots.txt

sed -i -e "s/.png/_-_${simple_name}.png/g" snapshots.txt
sed -i -e "s/snapshot /snapshot ${type}/g" snapshots.txt

echo "IGV.Bounds=0,0,1920,1080" > prefs.properties
echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties

xvfb-run --auto-servernum -s "-screen 0 1920x1080x24" java -Xmx2000m --module-path=/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties