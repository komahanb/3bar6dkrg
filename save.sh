read -p "Enter kprob value:" i
cp Opt.his kprob$i/
cp HISTG* kprob$i/
cp screen kprob$i/
cp fort.* kprob$i/
echo "Success saving the results into kprob"$i
