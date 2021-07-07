# Get bigwig:
mkdir -p input_pgt
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150701/suppl/*rev.bw" -P input_pgt

# Make the ini_file
gtf_file=$(ls input | grep gtf)
ini_file=input_pgt/Fgfr2.ini
echo "[x-axis]
where = top
" > ${ini_file}
for f in input_pgt/*MES*rev.bw input_pgt/*ECT*rev.bw; do
  title=$(basename $f .bw | awk -F "_" '{print $3}')
  if [[ $title = *"MES"* ]]; then
    color="#17B7FF"
  else
    color="#FC7C76"
  fi
  echo "[$title]
file = $(basename $f)
title = $title
min_value = 0
number_of_bins = 2000
show_data_range = false
color = $color
" >> ${ini_file}
done
echo "[spacer]
height = 0.5

[genes]
file = ../input/${gtf_file}
title = all exons merged
color_utr = black
merge_transcripts = true
prefered_name = gene_name
arrowhead_included = true
file_type = gtf

[spacer]
height = 0.5

[genes]
file = ../input/${gtf_file}
title = all transcripts
color_utr = black
arrowhead_included = true
height = 8
file_type = gtf
" >> ${ini_file}

# Plot Fig2C:
region=$(zcat input/${gtf_file} | grep "gene_name \"Fgfr2\"" | awk 'BEGIN{start=1000000000;end=1}{chr=$1;if($4<start){start=$4};if($5>end){end=$5}}END{printf("%s:%d-%d\n",chr,start,end + (end - start) / 10)}')
# chr7:130162451-130276644
pgt --tracks ${ini_file} --region ${region} -o figure/Fig2C.pdf --trackLabelFraction 0.1
