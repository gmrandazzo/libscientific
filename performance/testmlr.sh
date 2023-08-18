function runstuff(){
    row=$1
    col=$2
    memusage ./mlr_tests ${row} ${col} &> tmp
    m=`grep "heap total:" tmp  | cut -d ":" -f 3 | cut -d "," -f 1`
    t=`grep "Total time" tmp  | cut -d ":" -f 2 | cut -d " " -f 2`
    echo ${row},${col},${m},${t}
}

echo "row,cols,memory(bytes),time(sec)"

for row in `seq 0 100 5000`; do
    for col in `seq 0 200 1000`; do
        if [ ${row} -gt "0" ] && [ ${col} -gt "0" ]; then
            runstuff ${row} ${col}
        fi
    done
done


for row in `seq 0 5000 100000`; do
    for col in `seq 0 200 1000`; do
        if [ ${row} -gt "0" ] && [ ${col} -gt "0" ]; then
            runstuff ${row} ${col}
        fi
    done
done
