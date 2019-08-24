# ATAC-seq

## 必要なソフトウェアのインストール

### Picard のインストール

　PicardはJava 1.8に依存しているため、正常に動作するにはJava 1.8がインストールされている必要がある。まず、インストールされているjavaのバージョンを確認する。 java -version を実行して "1.8."で始まるバージョン名が表示されれば問題ない。

```
$ java -version
java version "1.8.0_45"
Java(TM) SE Runtime Environment (build 1.8.0_45-b14)
Java HotSpot(TM) 64-Bit Server VM (build 25.45-b02, mixed mode)
```
　Java 1.8がインストールされていない場合は、https://www.oracle.com/technetwork/java/javase/downloads/index.html からインストーラーをダウンロードしてインストールする。

　次に、以下のコマンドでpicardをインストールする。

```
$ brew install picard-tools
```

### ATAC-seqデータのダウンロード
　ヒトのGM12878細胞に対するATAC-seq実験を行ったデータを用いる（参考文献６）。通常、ATAC-seqではペアエンドシーケンシングが行われる。そのため、１つのサンプルについて、リードペアのRead 1およびRead 2に対応する２つのFASTQファイルをダウンロードする。Read 1とRead 2は\*\_1.fastq、\_\2.fastqというようなファイル名になっている場合が多い。
　まず、以下のコマンドでATAC-seq解析用のディレクトリを作成する。

```
$ mkdir ~/atacseq
```

　次に、ATAC-seqのFASTQファイルを格納するディレクトリを作成する。
```
$ mkdir ~/atacseq/fastq
```

　以下のコマンドで、ATAC-seqのFASTQファイルをダウンロードする。
```
$ cd ~/atacseq/fastq
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_1.fastq.gz
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz
```
　以下のファイルが~/atacseq/fastqにダウンロードされたことを確認する。
```
SRR891269_1.fastq.gz
SRR891269_2.fastq.gz
```

## ATAC-seq解析
### FASTQファイルのQC・前処理
#### FastQCによるFASTQファイルのQC
　以下のコマンドで、Read 1とRead 2のFASTQファイルそれぞれに対してFastQCを実行する。
```
$ cd ~/atacseq/
mkdir fastqc
fastqc -o fastqc data/fastq_atacseq/SRR891269_1.fastq.gz
fastqc -o fastqc data/fastq_atacseq/SRR891269_2.fastq.gz
```

　以下のファイルができたことを確認する。
```
SRR891269_1_fastqc.html
SRR891269_1_fastqc.zip
SRR891269_2_fastqc.html
SRR891269_2_fastqc.zip
```

　コマンドが終了したら、以下のコマンドでレポートを確認する。

```
$ cd ~/atacseq/
open fastqc/SRR891269_1_fastqc.html
open fastqc/SRR891269_2_fastqc.html
```

#### fastpによるリードトリミング
　以下のコマンドで、fastpによってFASTQのリードのトリミングを行う。ペアエンドリードの場合は、入力（--in1および--in2）と出力（--out1および--out2）についてそれぞれRead 1とRead 2に対応するファイル名を指定する。
```
$ cd ~/atacseq/
$ mkdir fastp
$ fastp --in1 data/fastq_atacseq/SRR891269_1.fastq.gz --in2 data/fastq_atacseq/SRR891269_2.fastq.gz \
	--out1 fastp/SRR891269_1.trim.fastq.gz --out2 fastp/SRR891269_2.trim.fastq.gz \
    --html fastp/SRR891269.fastp.html
```

　以下の出力ファイルがきちんと出力されているかを確認する。ペアエンドリードの場合、レポート（.htmlおよび.json）のファイルは１つになる。
```
fastp/SRR891269.fastp.html
fastp/SRR891269_1.trim.fastq.gz
fastp/SRR891269_2.trim.fastq.gz
```

　fastpのレポートを確認する。
```
$ cd ~/atacseq/
$ open fastp/SRR891269.fastp.html
```

### Bowtie2 によるペアエンドリードのマッピング
#### Macのコア数の確認
　以下のコマンドで、Macのコア数を確認する。出力結果はコア数が２の場合の例を示しており、実際には使用するマシンに依存して変化する。
```
$ system_profiler SPHardwareDataType | grep Cores
      Total Number of Cores: 2
```

#### Bowtie 2の実行
　Bowtie2によって、ATAC-seqのリードをヒトリファレンスゲノム配列にマッピングする。
　以下のコマンドで、Bowtie 2を実行する。-1 と -2 でRead 1、Read 2 のFASTQファイルを指定する。

```
$ cd ~/atacseq/
bowtie2 -p 2 --no-mixed --no-discordant -X 2000 \
-x data/external/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index \
-1 fastp/SRR891269_1.trim.fastq.gz -2 fastp/SRR891269_2.trim.fastq.gz > bowtie2/SRR891269.trim.sam
-X 2000：ペアエンドリードのマッピングの際に、Read 1の5’端とRead 2 の5’端の間の距離が2,000塩基以下になるマッピングだけを探索するように指定する。
--no-mixed：ペアエンドリードをペアとしてマッピングできなかった際に、シングルエンドリードとしてマップさせるようにしない。
--no-discordant：
```

　以下のコマンドで、SAMファイルからマップされなかったリードを除き、BAMに変換し、BAMをソートする。さらにBAMのインデックスを作成する。
```
$ cd ~/atacseq/
$ samtools view -bhS -F 0x4 bowtie2/SRR891269.trim.sam | \
samtools sort -T bowtie2/SRR891269.trim - \
> bowtie2/SRR891269.trim.bam
$ samtools index bowtie2/SRR891269.trim.bam
```

　以下のコマンドで、適切なペアとしてマップされたリード（"read mapped in proper pair"）のみを抽出する。
```
$ cd ~/atacseq/
$ samtools view -f 0x2 -bh bowtie2/SRR891269.trim.bam \
> bowtie2/SRR891269.trim.proper_pairs.bam
$ samtools index bowtie2/SRR891269.trim.proper_pairs.bam
$ picard MarkDuplicates \
    I=bowtie2/SRR891269.trim.proper_pairs.bam \
    O=bowtie2/SRR891269.trim.proper_pairs.rmdup.bam \
    M=bowtie2/SRR891269.trim.proper_pairs.rmdup.bam.log.txt \
    REMOVE_DUPLICATES=true
$ samtools index bowtie2/SRR891269.trim.proper_pairs.rmdup.bam
```

### Picardによるフラグメント長（インサートサイズ）の分布の可視化
　ATAC-seqは通常ペアエンドでシーケンシングされる。まず、ATAC-seqの実験がうまくいったかの検証のため、DNA断片の長さの分布を調べる。

```
$ cd ~/atacseq/
$ mkdir picard
$ picard CollectInsertSizeMetrics \
INPUT=bowtie2/SRR891269.trim.proper_pairs.rmdup.bam \
OUTPUT=picard/insert_size_metrics.txt \
HISTOGRAM_FILE=picard/hist.pdf \
MINIMUM_PCT=0
```

### MACS2によるピーク検出
　ChIP-seqの場合と異なり、ATAC-seqではインプットDNAなどのコントロール実験が行われないが多い。以下では、コントロールデータなしでATAC-seqからピーク検出を行う。
　以下のコマンドで、MACS2によるピーク検出を行う。
```
$ cd ~/atacseq
$ 	source ~/tools/MACS2/bin/activate   # MACS2がインストールされた環境へ切り替える
$ macs2 callpeak -t bowtie2/SRR891269.trim.proper_pairs.rmdup.bam \
-f BAM --nomodel --shift -50 --extsize 100 -g hs -n SRR891269_atacseq \
-B --outdir macs2
$ deactivate    # 元の環境へ切り替える
```

> TIPS：ATAC-seqデータにMACS2を使うときのコツ
? 　ATAC-seqではリードの5’端の位置（の4〜5塩基下流）がTn5の挿入部位に相当する。Tn5の挿入部位をピークとしてMACS2に検出させたい場合は、以下のようにリード長に依存したパラメータを設定する。この設定の仕方は、あたかもリードの5’端の位置を中心とした 2 x (リード長)の幅のリードがあるようにMACS2に認識させることに相当する。
> --shift (-1) x (リード長)
> --extsize 2 x (リード長)
>　上の例では、FastQCやFastpの結果からリード数が概ね50 bpであることがわかっているため、”--shift -50 --extsize 100”とした。


　以下のファイルが出力されることを確認する。
```
macs2/SRR891269_atacseq_peaks.narrowPeak  # ピーク領域を示すbedファイル
macs2/SRR891269_atacseq_peaks.xls    # ピーク領域と、関連するリード数等の詳細データ
macs2/SRR891269_atacseq_summits.bed    # ピーク領域での頂上部分を表すbedファイル
macs2/SRR891269_atacseq_treat_pileup.bdg    # リードの集積度を示すbedGraphファイル
macs2/SRR891269_atacseq_control_lambda.bdg    # バックグラウンドとして使用した、ポワソン分布におけるlambda値を表すbedGraphファイル
```

　以下のコマンドで、ピーク数を確認する。
```
$ cd ~/atacseq
$ wc -l macs2/SRR891269_atacseq_peaks.narrowPeak
   51129 macs2/SRR891269_atacseq_peaks.narrowPeak
	 ```

　以下のコマンドで、染色体ごとのピーク数を確認する。
```
$ cd ~/atacseq
$ cut -f 1 macs2/SRR891269_atacseq_peaks.narrowPeak | uniq -c | sort -r
4775 chr1
4500 chr2
3633 chr3
3405 chr6
3243 chr5
3147 chr4
(以下省略)
```

### IGV によるピークの確認
　ChIP-seqの時と同様に、ポジティブコントロールの領域においてMACS2による検出されたピークが存在するかをIGVで確認する。
　以下のコマンドで、IGVを起動する。
```
$ igv
```
