# 162个宏基因组计算life-history trait values

## 1.计算平均基因组大小（ags)

批量更改后缀，将fa改成fna/faa

```
for file in *.fa; do mv "$file" "${file%.fa}.fna"; done
for file in *.fa; do mv "$file" "${file%.fa}.faa"; done
```

批量运行ags计算：

#将下列脚本存储为run_samples.sh

```
#!/bin/bash

samples=("S9" "S72" "S74" "S95" "S61" "S161")  # 更新为你的样本名列表

for sample in ${samples[@]}; do
  fna_file="//Media/E/fanqilian/fna/${sample}.fna"
  faa_file="//Media/E/fanqilian/faa/${sample}.faa"
  out_file="${sample}.out"

  nohup /home/qilianfan/mg/AGS-and-ACN-tools/run_ags.sh \
    $fna_file \
    $faa_file \
    /home/qilianfan/mg/AGS-and-ACN-tools/out/${sample} \
    --sample_name $sample \
    --verbose t \
    --min_length 120\
    --overwrite t \
    --nslots 4 \
    --train_file_name complete \
    --save_complementary_data t >> $out_file &
done
```

然后运行

```
bash run_samples.sh
```

```
find . -name "*_ags.tsv" -exec cp {} /home/qilianfan/mg/AGS-and-ACN-tools/out/tsv/ \;
```

```
import os
import pandas as pd

# 定义需要遍历的文件夹路径
folder_path = 'H:/样本的faa,fna,gff文件/ags/'

# 创建一个空的DataFrame来存储所有文件的数据
merged_data = pd.DataFrame()

# 遍历文件夹下的所有文件
for file_name in os.listdir(folder_path):
    # 检查文件扩展名是否为tsv
    if file_name.endswith('.tsv'):
        # 构建文件的完整路径
        file_path = os.path.join(folder_path, file_name)
        
        # 读取tsv文件数据
        data = pd.read_csv(file_path, sep='\t')
        
        # 将数据添加到合并的DataFrame中
        merged_data = merged_data.append(data)

# 去除重复行
merged_data.drop_duplicates(inplace=True)

# 将合并后的数据写入新的文件
output_file = 'H:/样本的faa,fna,gff文件/ags_output.tsv'  # 指定输出文件的路径和名称
merged_data.to_csv(output_file, sep='\t', index=False)

print("合并数据已写入文件:", output_file)
```



## 2.计算16S平均拷贝数(注意所有输入数据在同一目录下)

```
#!/bin/bash

samples=("S141" "S161")  # 更新为你的样本名列表
for sample in ${samples[@]}; do
  fna_file="/home/qilianfan/mg/AGS-and-ACN-tools/out/${sample}/${sample}_FBL.fna"
  tsv_file="/home/qilianfan/mg/AGS-and-ACN-tools/out/${sample}/${sample}_ags.tsv"
  out_file="${sample}.out"

  nohup /home/qilianfan/mg/AGS-and-ACN-tools/run_acn.sh \
    $fna_file \
    $tsv_file \
    /home/qilianfan/mg/AGS-and-ACN-tools/acn_out/${sample} \
    --sample_name $sample \
    --verbose t \
    --overwrite t \
    --nslots 4 \
    --save_complementary_data t >> $out_file &
done
```

```
#将每个样品文件夹下抽取.tsv文件，然后合并所有文件并去除重复行
import os
import pandas as pd

# 定义A文件夹的路径
folder_path = 'H:/样本的faa,fna,gff文件/acn/'

# 创建一个空的DataFrame来存储所有文件的数据
merged_data = pd.DataFrame()

# 遍历A文件夹下的每个子文件夹
for folder_name in os.listdir(folder_path):
    # 构建子文件夹的完整路径
    folder = os.path.join(folder_path, folder_name)
    
    # 检查路径是否为文件夹
    if os.path.isdir(folder):
        # 遍历子文件夹下的所有文件
        for file_name in os.listdir(folder):
            # 检查文件扩展名是否为tsv
            if file_name.endswith('.tsv'):
                # 构建文件的完整路径
                file_path = os.path.join(folder, file_name)
                
                # 读取tsv文件数据
                data = pd.read_csv(file_path, sep='\t')
                
                # 将数据添加到合并的DataFrame中
                merged_data = merged_data.append(data)

# 去除重复行
merged_data.drop_duplicates(inplace=True)

# 将合并后的数据写入新的文件
output_file = 'H:/样本的faa,fna,gff文件/acn_output.tsv'  # 指定输出文件的路径和名称
merged_data.to_csv(output_file, sep='\t', index=False)

print("合并数据已写入文件:", output_file)
```



## 3. 利用growpred计算Maximum growth rate (h-1）和Condon usage bias

首先利用growthpred计算Minimum growth rate然后取倒数

Growthsnake: Snakemake pipeline to run growthpred on large amount of genomes.

https://gitlab.univ-nantes.fr/combi-ls2n/growthsnake

```
git clone https://gitlab.univ-nantes.fr/combi-ls2n/growthsnake
cd growthsnake
tar -zxvf growthpred.tar.gz
```

因为官网是网页版，所以从growthsnake里得到growthpred的linux软件包，**Growthsnake跑不出来，所以直接使用Growthpred**

```
#!/bin/bash

mkdir -p test_output
./growthpred-v1.08.py -d ./ -f ecoli_ribosomal_genes.txt -g ecoli_complete_genome.txt -s -S -c 0 -o test_output/test_ecoli -t -m
./growthpred-v1.08.py -d ./ -g ecoli_complete_genome.txt -s -S -c 0 -o test_ecoli -t -m -b
./growthpred-v1.08.py -f PTR0Bin1Ribosomal.fasta -g PTR0Bin1nucgenes.fna -s -S -c 0 -o test_PTR0Bin1 -t -m
./growthpred-v1.08.py -f PTR1BBin29ribo.fasta -g PTR1BBin29nucgenes.fna -s -S -c 0 -o test_PTR1Bin29 -t -m

./growthpred-v1.08.py -d ./ -g TARA_023_SRF_redu2M.fna -s -S -c 0 -o test/test -t -m -b
```

运行前先对growthpred-v1.08.py脚本进行修改，重置temp文件夹路径。此外，可能环境中需要安装一下R

```
GROWTHPRED_DIR = os.path.dirname(os.path.abspath(__file__))
GROWTHPRED_BIN = GROWTHPRED_DIR + "/"
GROWTHPRED_LIBEX = os.path.join(GROWTHPRED_DIR, "Programs/")
GROWTHPRED_SHARE = os.path.join(GROWTHPRED_DIR, "shared/")
GROWTHPRED_TMP = "/lustre/home/acct-ioozy/ioozy-user3/work/growthsnake/growthpred/temp/"
GROWTHPRED_R = GROWTHPRED_BIN
```

在超算运行：

```
#!/bin/bash
#SBATCH -J growthpred
#SBATCH -p cpu
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=hahahehe@sjtu.edu.cn

for i in $(seq -f "S%g" 1 10); do
    /lustre/home/acct-ioozy/ioozy-user3/work/growthsnake/growthpred/growthpred-v1.08.py  -d /lustre/home/acct-ioozy/ioozy-user3/work/fna/ -g ${i}.fna -s -S -c 0 -o ${i} -t -m -b -r
done

```

### 3.1 计算Maximum growth rate (h-1）

跑完growthpred后，得到.results文件，对该文件进行以下处理：

```
import os
import re

# 路径到包含 results 文件的文件夹
directory = "G:/MG/样本的faa,fna,gff文件/MGT/"

# 创建一个空的列表来保存结果
results = []

# 对文件夹中的每个文件进行迭代
for filename in os.listdir(directory):
    if filename.endswith(".results"):
        with open(os.path.join(directory, filename), 'r') as file:
            for line in file:
                # 寻找包含预测的最小生成时间的行
                if line.startswith("Predicted minimum generation time"):
                    # 使用正则表达式提取小时数
                    time = re.search(r'(\d+\.\d+)', line).group()
                    # 添加文件名（不带扩展名）和生成时间到结果列表
                    results.append([filename.rstrip('.results'), time])

# 将结果写入新的 CSV 文件
with open('G:/MG/样本的faa,fna,gff文件/MGT/generation_times.csv', 'w') as file:
    for result in results:
        file.write(f"{result[0]},{result[1]}\n")
```

然后对得到的minimum generation time值求倒数就得到Maximum growth rate (h-1）

### 3.2 计算Condon usage bias

跑完growthpred后得到.cub文件，对该文件进行下列处理：

```
import os
import pandas as pd
import numpy as np

# 路径到包含 cub 文件的文件夹
directory = "G:/MG/样本的faa,fna,gff文件/ENC/"

# 创建一个空的 DataFrame 来保存结果
results = pd.DataFrame(columns=['filename', 'average_ENCp'])

# 对文件夹中的每个文件进行迭代
for filename in os.listdir(directory):
    if filename.endswith(".cub"):
        # 读取文件
        data = pd.read_csv(os.path.join(directory, filename), sep='\t')
        # 计算 ENCp 列的平均值，忽略 NA 值
        average_ENCp = data['ENCp'].replace('NA', np.nan).astype(float).mean()
        # 添加文件名（不带扩展名）和平均 ENCp 到结果 DataFrame
        results = results.append({'filename': filename.rstrip('.cub'), 'average_ENCp': average_ENCp}, ignore_index=True)

# 将结果写入新的 CSV 文件
results.to_csv('G:/MG/样本的faa,fna,gff文件/ENC/average_ENCp.csv', index=False)
```

求个每个样本所有核糖体基因的平均ENC‘后求倒数得到所要的Condon usage bias

## 4. 求GC含量和GC含量方差

```
#对单个文件
from Bio import SeqIO
import numpy as np

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

# 读取FASTA文件
sequences = list(SeqIO.parse("G:/MG/样本的faa,fna,gff文件/fna/S1.fa", "fasta"))

# 计算每个序列的GC含量
gc_contents = [gc_content(str(seq_record.seq)) for seq_record in sequences]

# 计算GC含量的平均值和方差
average_gc_content = np.mean(gc_contents)
variance_gc_content = np.var(gc_contents)

print(f'Average GC content: {average_gc_content}')
print(f'Variance of GC content: {variance_gc_content}')

```

```
#批量求取
import os
from Bio import SeqIO
import numpy as np

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

# 指定你的文件夹路径
folder_path = "G:/MG/样本的faa,fna,gff文件/fna/16/"

# 创建一个空的列表来存储结果
results = []

# 遍历文件夹下的所有.fa文件
for filename in os.listdir(folder_path):
    if filename.endswith('.fa'):
        sample_name = filename[:-3]  # remove .fa from filename for the sample name
        sequences = list(SeqIO.parse(folder_path + filename, "fasta"))
        gc_contents = [gc_content(str(seq_record.seq)) for seq_record in sequences]
        average_gc_content = np.mean(gc_contents)
        variance_gc_content = np.var(gc_contents)
        results.append((sample_name, average_gc_content, variance_gc_content))

# 输出结果到文件
with open("G:/MG/样本的faa,fna,gff文件/fna/16/GC_content_results.txt", "w") as output_file:
    output_file.write("Sample Name\tAverage GC Content\tVariance of GC Content\n")
    for result in results:
        output_file.write(f'{result[0]}\t{result[1]}\t{result[2]}\n')

```

