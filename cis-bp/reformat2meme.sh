ls pwms/ | grep txt | perl -e 'while(<>){chomp;$m=substr($_,0,-4);system("python reformat2meme.py -i ./pwms/$m.txt -m $m -o ./pwms/$m.meme");}'
