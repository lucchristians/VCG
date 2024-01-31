num=1
path=$1
name=$2
rm ${path}/${name}.txt
touch ${path}/${name}.txt
#echo "" > ${name}.txt
until [ ! -f ${path}/${name}${num}.txt ]
do
	cd ${path}
	cat ${name}.txt ${name}${num}.txt >> ${name}_temp.txt
	mv ${name}_temp.txt ${name}.txt
	((num++))
	cd ..
done	
