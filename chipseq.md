# ChIP-seq

## ソフトウェアのインストール

### fastpのインストール
　fastpはリードのトリミングを行うツールである。
　以下のコマンドでfastpをインストールする。

```
$ git clone https://github.com/OpenGene/fastp.git
$ cd fastp
$ make
$ sudo make install
$ fastp --help
```

### Bowtie 2 のインストール
#### Bowtie 2プログラムのインストール

```
$ conda install -c bioconda bowtie2
```

#### Pre-built indexの取得
　マッピングのためのデータベースPre-built indexを取得する。Bowtie 2 のページ（http://bowtie-bio.sourceforge.net/bowtie2/index.shtml ）へアクセスする。右側の”Indexes”という箇所から、対象の生物種に合ったPre-built indexをダウンロードする。

　本項で使用する種はマウスなので、M. Musculus (mm10) のリンク先アドレス（ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip）をコピーし、以下のコマンドでマウスゲノム (mm10)のインデックスファイルをダウンロードし、解凍する。

```
$ mkdir ~/bowtie2_index    # Pre-built indexを入れるためのディレクトリを作成する
$ cd ~/bowtie2_index    # # Pre-built indexを入れるためのディレクトリに移動する
$ wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
$ unzip mm10.zip
```
　同様に、以下のコマンドで、ヒトリファレンスゲノム配列（hg38）用のpre-built indexをダウンロードし、さらにダウンロードされた圧縮ファイルを解凍する。

```
$ cd ~/bowtie2_index
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
$ tar xvzf bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
```

　なお、IlluminaのiGenomesのページ (https://support.illumina.com/sequencing/sequencing_software/igenome.html) においても様々な生物種のゲノムに対するpre-built indexが配布されている。

### MACS2のインストール
　現行のMACS2は標準ではpython 2 のみに対応している。一方で、python 2 は2020年1月にサポートが終了する（参考文献７）。今後の継続性を考慮し、ここではpython 3 でのインストール方法を紹介する。なお、2019年2月現在では、 `pip` や `conda` といった はpython 2のみで正常にインストールできる。
　まず、MACS2をインストールする環境を作る。

```
$ mkdir ~/tools
$ cd ~/tools
$ python3 -m venv MACS2/
$ source MACS2/bin/activate
```

　以下のコマンドでMACS2のpython 3版をインストールする。

```
$ pip install --upgrade pip
$ pip install numpy
$ pip install cython
$ git clone https://github.com/taoliu/MACS.git
$ cd MACS
$ git checkout remotes/origin/macs2python3
$ python setup_w_cython.py install
```

　以下のコマンドでMACS2がインストールされたことを確認する。

```
$ macs2 --help  # ヘルプメッセージが出力されればOK
```

　以下のコマンドで、MACS2をインストールした環境から離脱する。

```
$ deactivate
```

　MACS2を使用する際には、あらかじめ以下のコマンドでMAC2をインストールした環境に変えればよい。
```
$ source ~/tools/MACS2/bin/activate
```

### samtools のインストール
```
$ brew install samtools
```

### HOMERのインストール
#### HOMERのインストール
```
$ conda install -c bioconda homer
```
#### HOMERで特定のゲノムを使えるようにする

　以下のコマンドでHOMERでマウスゲノム（hg38）を使えるようにする。
```
$ perl /anaconda3/share/homer-*/configureHomer.pl -install hg38
```

　以下のコマンドでHOMERでマウスゲノム（mm10）を使えるようにする。

```
$ perl /anaconda3/share/homer-*/configureHomer.pl -install mm10
```

### deepToolsのインストール
```
$ conda install -c bioconda deeptools samtools=1.9=h8ee4bcc_1 openssl=1.0
$ deeptools --version
```

### RおよびRStudioのインストール
#### Rのインストール
```
$ brew install r
```

#### RStudioのインストール
```
$ brew cask install rstudio
```

### IGVのインストール
```
$ brew install igv
```

### bedtools のインストール

```
$ brew install bedtools
```

### ChIPpeakAnno のインストール

　まず、以下のコマンドで、RStudioを起動する。

```
$ open -a RStudio
```

次に、以下のRのコマンドで、chipPeakAnnoをインストールする。なお、"Update all/some/none? [a/s/n]:" と表示されたら "a" と入力してエンターを押すようにする。

```
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("ChIPpeakAnno")
```

## 使用するデータの取得
### 遺伝子モデルのダウンロード

- GENCODE https://www.gencodegenes.org
    - ヒトの場合は https://www.gencodegenes.org/human/
    - マウスの場合は https://www.gencodegenes.org/mouse/

上記URLにアクセスし、Contentが"Comprehensive gene annotation"、Regionsが"CHR"である行のGTFのダウンロードリンクをコピーする。


　以下のコマンドで、GENCODEのマウスのリファレンス遺伝子モデルをダウンロードする。
```
$ mkdir ~/gencode
$ cd ~/gencode
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz # このURLはコピーしたダウンロードリンクに応じて変わる
$ gzip -d gencode.vM20.annotation.gtf.gz
$ ls gencode.vM20.annotation.gtf　# gencode.vM20.annotation.gtfができたことを確認する。
```

### 使用するChIP-seqのFASTQファイルのダウンロード

　以下のコマンドで、ChIP-seq解析用のディレクトリを作成し、さらにその下にFASTQファイルを保存するディレクトリを作成する。

```
$ mkdir ~/chipseq    # ChIP-seq解析用のディレクトリを作成する
$ cd ~/chipseq    # ChIP-seq解析用のディレクトリへ移動する
$ mkdir fastq    # FASTQファイルを入れるディレクトリを作成する
```

　マウスで行われたChIP-seq実験のFASTQファイルをダウンロードする。
　ここでは、EMBL-EBIが運営するENA (European Nucleotide Archive) からFASTQファイルをダウンロードする。

1）マウスAT-3細胞（IFN-γ添加時）におけるBRD4 ChIP-seq
　SRR5208824.fastq.gzをダウンロードする。

```
$ cd ~/chipseq/fastq
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR520/004/SRR5208824/SRR5208824.fastq.gz
```

2）マウスAT-3細胞（IFN-γ添加時）におけるIRF1 ChIP-seq
　SRR5208828.fastq.gzをダウンロードする。

```
$ cd ~/chipseq/fastq
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR520/008/SRR5208828/SRR5208828.fastq.gz
```

３）マウスAT-3細胞におけるInput DNA ChIP-seq
Input DNA実験はChIP-seqの実験のネガティブコントロールとして用いらる。ChIP-seqのコントロール実験のサンプルは、MACS2によるピーク検出の際に非特異的なピークを除去するために用いられる。
　SRR5208838.fastq.gzをダウンロードする。

```
$ cd ~/chipseq/fastq
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR520/008/SRR5208838/SRR5208838.fastq.gz
```


## ChIP-seq解析
### FASTQファイルの名称変更

　SRR5208824.fastq.gzをコピーし、BRD4_ChIP_IFNy.R1.fastq.gzに名前を変更する

```
$ cd ~/chipseq/fastq
$ cp SRR5208824.fastq.gz BRD4_ChIP_IFNy.R1.fastq.gz
```

　以下、他の２つのファイルについても同様に実行する。
```
$ cp SRR5208828.fastq.gz IRF1_ChIP_IFNy.R1.fastq.gz
$ cp SRR5208838.fastq.gz Input_DNA.R1.fastq.gz
```

### FASTQファイルのQC・前処理
#### FastQCによるリードのQC
　まずFastQCの結果を入れるディレクトリを作成する。

```
$ cd ~/chipseq
$ mkdir fastqc
```

　以下のコマンドで、FastQCを実行する。
```
$ cd ~/chipseq
$ fastqc -o fastqc fastq/BRD4_ChIP_IFNy.R1.fastq.gz
```

　同様に、他のFASTQファイルに対してもFastQCを実行する。
```
$ cd ~/chipseq
$ fastqc -o fastqc fastq/IRF1_ChIP_IFNy.R1.fastq.gz
$ fastqc -o fastqc fastq/Input_DNA.R1.fastq.gz
```

　fastqcディレクトリの中に、以下のファイルができていることを確認する。
```
BRD4_ChIP_IFNy.R1_fastqc.html
BRD4_ChIP_IFNy.R1_fastqc.zip
IRF1_ChIP_IFNy.R1_fastqc.html
IRF1_ChIP_IFNy.R1_fastqc.zip
Input_DNA.R1_fastqc.html
Input_DNA.R1_fastqc.zip
```

#### fastp によるFASTQファイルのリードトリミング
　まずfastpの結果を入れるディレクトリを作成する。

```
$ cd ~/chipseq
$ mkdir fastp    # fastpという名前のディレクトリを作る
```

　以下のコマンドで、BRD4_ChIP_IFNy.R1.fastq.gzに対してfastpを実行する。ここではリードトリミング後のFASTQファイルをBRD4_ChIP_IFNy.R1.trim.fastq.gzとして出力させる。

```
$ cd ~/chipseq
$ fastp -i fastq/BRD4_ChIP_IFNy.R1.fastq.gz -o fastp/BRD4_ChIP_IFNy.R1.trim.fastq.gz --html fastp/BRD4_ChIP_IFNy.fastp.html
```

　同様にIRF1_ChIP_IFNy.R1.fastq.gz、Input_DNA.R1.fastq.gzに対してfastpを実行する。
```
$ cd ~/chipseq
$ fastp -i fastq/IRF1_ChIP_IFNy.R1.fastq.gz -o fastp/IRF1_ChIP_IFNy.R1.trim.fastq.gz --html fastp/IRF1_ChIP_IFNy.fastp.html
$ fastp -i fastq/Input_DNA.R1.fastq.gz -o fastp/Input_DNA.R1.trim.fastq.gz --html fastp/Input_DNA.fastp.html
```

　fastpディレクトリ内に以下のファイルができてることを確認する。

```
BRD4_ChIP_IFNy.fastp.html
BRD4_ChIP_IFNy.R1.trim.fastq.gz
IRF1_ChIP_IFNy.fastp.html
IRF1_ChIP_IFNy.R1.trim.fastq.gz
Input_DNA.fastp.html
Input_DNA.R1.trim.fastq.gz
```

　以下のコマンドで、fastpのレポートを確認できる。
```
$ cd ~/chipseq
$ open fastp/BRD4_ChIP_IFNy.fastp.html
```


> Tips：FASTQのリードの3’端の塩基をすべてのリードから除く操作が必要な場合がある
　Illumina社製のDNAシーケンサーの最後のサイクル（最後の塩基）では塩基読み取り精度が低くなることが知られている。そのため、シーケンサーを動かす実験者やシーケンシングセンター、受託企業の方針によっては、目的のリード長のサイクルに１塩基分サイクルを追加し、FASTQファイルの前処理の段階で3’端の１塩基を削る場合がある。例えば、100塩基長のリードがほしいときには、101サイクル回してリード長が101塩基のリードからなるFASTQファイルを取得し、データ解析の段階で3’端の１塩基を削って100の塩基のリード長のFASTQファイルを出力させるといった具合である。
>　このような場合には、fastpなどのリードトリミングツールによって最後の１塩基を削る操作が必要となる。fastpではすべてのリードが同じ長さの場合、`--trim_tail1=1`というオプションを使用することで、すべてのリードの3'端の1塩基を削ることができる（ペアエンドリードの場合は`--trim_tail1=1 --trim_tail2=1`）
>　一方、すでにリードトリミング後のFASTQファイルである場合、 `--max_len1 N` (Nは正の整数) というオプションを使用することで、Read 1でN塩基を越えるリードがあったらN塩基になるまで3'端から塩基を削るという処理を行うことができる（ペアエンドリードの場合は `--max_len1 N --max_len2 N`）。

### リードのゲノムへのマッピング
#### Macのコア数の確認
　Bowtie 2によるリファレンスゲノム配列へのリードのマッピングは計算が重いため、計算終了までに数時間から半日かかる場合がある。そこで、複数のコアに計算を分散させることで計算の高速化を図りたい。そのためにまず、Macのコア数を確認する。
　以下のコマンドで、Macのコア数を確認する。以下の例では、Total Number of Cores: 2と表示されている。最近のMacで搭載されているIntelのCoreシリーズではハイパースレッディング・テクノロジー (Hyper-Threading Technology)が用いられているため、１つのコアで２スレッドが動作するため、同時に動作可能なスレッド数は4となる。

```
$ system_profiler SPHardwareDataType | grep Cores
      Total Number of Cores: 2
```

#### Bowtie 2によるリードのゲノムへのマッピング
　次に、Bowtie 2によってリードをリファレンスゲノム配列へマッピングする。
　まず、Bowtie 2の計算結果を入れるディレクトリを作成する。

```
$ cd ~/chipseq
$ mkdir bowtie2
```

　以下のコマンドで、Bowtie 2により、BRD4_ChIP_IFNy.R1.trim.fastq.gzのリードをマウスリファレンスゲノム（mm10）へマッピングする。

```
$ cd ~/chipseq
$ bowtie2 -p 2 -x data/external/bowtie2_index/mm10 \
    -U BRD4_ChIP_IFNy.R1.trim.fastq.gz > bowtie2/BRD4_ChIP_IFNy.trim.sam
```

> -p：使用するコア数を指定する。
> -x：リファレンスゲノム配列のPre-built indexの接頭辞を指定する。
> -U：FASTQファイルを指定する。
> なお、Bowtie 2 はデフォルトでは最もスコアが高いアラインメントを１つだけ出力する。もし最もスコアが高いアラインメントが複数見つかった場合は、その中からランダムに一つ選ぶ。

　計算が終了する際に、以下のようなどのくらいのリードがマッピングされたかの割合が出力される。一般に、overall alignment rateが極端に低い場合は、実験がうまくいっていなかったり、リファレンスゲノム配列の選択が不適切であるなど何らかの異常の可能性がある。

```
19709457 reads; of these:
  19709457 (100.00%) were unpaired; of these:
    1163507 (5.90%) aligned 0 times
    13206029 (67.00%) aligned exactly 1 time
    5339921 (27.09%) aligned >1 times
94.10% overall alignment rate
```

　bowtie2フォルダの中に、BRD4_ChIP_IFNy.trim.samが出力されていることを確認する。これはSAM形式のファイルと呼ばれる。

　同様に、他のFASTQファイルに対してもBowtie 2を実行する。

```
$ cd ~/chipseq
$ bowtie2 -p 2 -x data/external/bowtie2_index/mm10 \
    -U IRF1_ChIP_IFNy.R1.trim.fastq.gz > bowtie2/IRF1_ChIP_IFNy.trim.sam
$ bowtie2 -p 2 -x data/external/bowtie2_index/mm10 \
    -U Input_DNA.R1.trim.fastq.gz > bowtie2/Input_DNA.trim.sam
```

　以下のファイルができたことを確認する。

```
bowtie2/IRF1_ChIP_IFNy.trim.sam
bowtie2/Input_DNA.trim.sam
```

　SAMファイルを直接扱うよりは、SAMファイルをバイナリ化して圧縮したBAMファイルへ変換した方が後々の解析で便利である。そこで、samtoolsを用いてBowtie 2から出力されたSAMファイルをBAMファイルへ変換する。
　以下のコマンドでは、(1) SAMをBAMに変換し、(2) SAMからユニークなリード（複数のゲノム領域にマップされたリード）を抽出し、(3) さらにBAMをソートしている。

```
$ cd ~/chipseq
$ samtools view -bhS -F 0x4 -q 42 bowtie2/BRD4_ChIP_IFNy.trim.sam | samtools sort -T bowtie2/BRD4_ChIP_IFNy.trim - > bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam
```

`-F 0x4`によってマップされなかったリードを除き、`-q 42` によってユニークなリードだけを抽出することができる。

　以下のファイルができたことを確認する。
```
bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam
```

　同様に、残りのSAMファイルについてもBAMへの変換を行う。
```
$ cd ~/chipseq
$ samtools view -bhS -F 0x4 -q 42 bowtie2/IRF1_ChIP_IFNy.trim.sam | samtools sort -T bowtie2/IRF1_ChIP_IFNy.trim - > bowtie2/Input_DNA.trim.uniq.bam
$ samtools view -bhS -F 0x4 -q 42 bowtie2/Input_DNA.trim.sam | samtools sort -T bowtie2/Input_DNA.trim - > bowtie2/Input_DNA.trim.uniq.bam
```

　BAMファイルを読み込む際にBAMのインデックス (拡張子が.bam.bai) が必要になる場合が多い。そこで、以下のコマンドで、BAMのインデックスを作成する。

```
$ cd ~/chipseq
$ samtools index bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam
$ samtools index bowtie2/IRF1_ChIP_IFNy.trim.uniq.bam
$ samtools index bowtie2/Input_DNA.trim.uniq.bam
```

　以下のファイルができたことを確認する。

```
bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam.bai
bowtie2/IRF1_ChIP_IFNy.trim.uniq.bam.bai
bowtie2/Input_DNA.trim.uniq.bam.bai
```

### MACS2によるピーク検出

```
$ cd ~/chipseq
$ mkdir macs2    # MACS2の出力結果を保存するディレクトリを作成する（必須ではない）
```

　以下のコマンドで、bowtie2/BRD4_ChIP_IFNy.trim.uniq.bamおよびbowtie2/IRF1_ChIP_IFNy.trim.uniq.bamに対してそれぞれMACS2を適用し、ピーク検出を行う。なお、ここでは、bowtie2/Input_DNA.trim.uniq.bamをChIP-seq実験のネガティブコントロールとして “-c” で指定している。

```
$ cd ~/chipseq
$ source ~/tools/MACS2/bin/activate   # MACS2がインストールされた環境へ切り替える
$ macs2 callpeak -t bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam \
-c bowtie2/Input_DNA.trim.uniq.bam -f BAM -g mm -n BRD4_ChIP_IFNy --outdir macs2 -B -q 0.01
$ macs2 callpeak -t bowtie2/IRF1_ChIP_IFNy.trim.uniq.bam \
-c bowtie2/Input_DNA.trim.uniq.bam -f BAM -g mm -n IRF1_ChIP_IFNy --outdir macs2 -B -q 0.01
$ deactivate    # 元の環境へ切り替える。
```

> -g：生物種を指定する。マウスだとmm、ヒトだとhsにする。
> -q：ピークを出力する際の「補正されたp値」（adjusted p-value）あるいはq値(q-value)の閾値を表す。デフォルトではq値の閾値は0.05である。p値ではなくq値を用いるのは、多重検定補正のためである。MACS2ではピークの”確からしさ”に関する統計的仮説検定をピーク候補の数だけ行う（多重検定）。このような場合、たとえ帰無仮説を棄却できない場合でも何度も検定を行えば、少なくとも１度は帰無仮説が棄却される割合が検定回数に従って増え、偽陽性の危険性が高まる。そのため、多重検定補正が必要となる。この場合、通常のp値ではなくp値に多重検定補正を施したq値(q-value)やFDR(False discovery rate)で閾値を設定する。MACS2では、BH法 (Benjamini-Hochberg method)を用いてFDRを計算している。多重検定補正については参考文献1を参照されたい。
> -B：bigBedファイルを出力させる
> -c：コントロールのデータを指定する。ChIP-seqのコントロール実験が実施されていないといった理由からコントロールのデータを用いない場合は使用しない。

　以下のファイルができていることを確認する。
```
BRD4_ChIP_IFNy_model.r  # ChIP DNAフラグメント長の推定結果
BRD4_ChIP_IFNy_peaks.narrowPeak  # ピーク領域を表すnarrowPeakファイル
BRD4_ChIP_IFNy_peaks.xls    # ピークの詳細情報を示すExcelファイル
BRD4_ChIP_IFNy_summits.bed    # ピーク領域の中で頂上となる部分 (summit) を表すBEDファイル
BRD4_ChIP_IFNy_treat_pileup.bdg   #
BRD4_ChIP_IFNy_control_lambda.bdg   #
```

　peaks_.narrowPeak はBED6+4 format形式のファイルである。  \*\_peaks.narrowPeak や \*\_summits.bed では、検出されたピークの情報が１行ずつ記載されている。

　例えば、以下のheadコマンドを使って、ファイルの最初の10行を表示することができる。
```
$ head macs2/BRD4_ChIP_IFNy_peaks.narrowPeak
chr1	4807514	4808176	BRD4_ChIP_IFNy_peak_1	117	.	8.68079	15.22217	11.75610	520
chr1	4857437	4857680	BRD4_ChIP_IFNy_peak_2	35	.	4.76915	6.42381	3.56659	44
chr1	4857758	4858397	BRD4_ChIP_IFNy_peak_3	216	.	12.20127	25.56099	21.63464	266
chr1	5018884	5019146	BRD4_ChIP_IFNy_peak_4	48	.	5.54587	7.84017	4.84834	123
chr1	5019310	5019671	BRD4_ChIP_IFNy_peak_5	60	.	6.07428	9.09495	6.02557	165
chr1	5022794	5023366	BRD4_ChIP_IFNy_peak_6	117	.	8.68079	15.22217	11.75610	450
chr1	5082960	5083202	BRD4_ChIP_IFNy_peak_7	80	.	5.95644	11.26737	8.04385	189
chr1	6214644	6215164	BRD4_ChIP_IFNy_peak_8	95	.	7.50730	12.85034	9.50346	161
chr1	7088391	7088703	BRD4_ChIP_IFNy_peak_9	73	.	6.72497	10.58313	7.39270	149
chr1	9747880	9748331	BRD4_ChIP_IFNy_peak_10	78	.	6.70945	11.01250	7.80631	189
```
　\*\_peaks.narrowPeakは

> １列目：染色体番号
> ２列目：ピークの５'端
> ３列目：ピークの３’端
> ４列目：ピーク名
> ５列目：ピークの-10*log10(qvalue)を整数に変換した値
> ６列目：ピークのストランド（ChIP-seqのピークはストランド情報は無いため.と表示）
> ７列目：バックグラウンドとのfold change
> ８列目：-log10(pvalue)
> ９列目：-log10(qvalue)
> 10列目：ピークの５'端からピークの頂上への相対的な位置

　また、\*\_summits.bedは検出されたピークの頂上部分の位置を示す。列の説明は \*\_peaks.narrowPeak の１〜５行目に相当する。
```
$  head macs2/BRD4_ChIP_IFNy_summits.bed
chr1	4808034	4808035	BRD4_ChIP_IFNy_peak_1	11.75610
chr1	4857481	4857482	BRD4_ChIP_IFNy_peak_2	3.56659
chr1	4858024	4858025	BRD4_ChIP_IFNy_peak_3	21.63464
chr1	5019007	5019008	BRD4_ChIP_IFNy_peak_4	4.84834
chr1	5019475	5019476	BRD4_ChIP_IFNy_peak_5	6.02557
chr1	5023244	5023245	BRD4_ChIP_IFNy_peak_6	11.75610
chr1	5083149	5083150	BRD4_ChIP_IFNy_peak_7	8.04385
chr1	6214805	6214806	BRD4_ChIP_IFNy_peak_8	9.50346
chr1	7088540	7088541	BRD4_ChIP_IFNy_peak_9	7.39270
chr1	9748069	9748070	BRD4_ChIP_IFNy_peak_10	7.80631
```



　次に、いくつのピークが検出されたかを数える。ピーク１つが１行で表せされるので、wc -l でファイルの行数を計算すれば、検出されたピーク数がわかる。

```
$ cd ~/chipseq
$ $ wc -l macs2/*_peaks.narrowPeak
    9348 macs2/BRD4_ChIP_IFNy_peaks.narrowPeak
     907 macs2/BRD4_ChIP_IFNy_peaks.overlapped_with_IRF1_ChIP_IFNy_peaks.narrowPeak
    3866 macs2/IRF1_ChIP_IFNy_peaks.narrowPeak
```

　次に、BRD4のピークとIRF1のピークがどのくらい重なるかを調べる。２つのピーク集合の間での重なりを調べるためにbedtoolsを使用する。
　以下のコマンドでは、-a で指定したピーク群（BRD4）のうち、-bで指定したピーク群（IRF1）と重なるものを抽出する。

```
$ cd ~/chipseq
$ bedtools intersect -u -a macs2/BRD4_ChIP_IFNy_peaks.narrowPeak -b macs2/IRF1_ChIP_IFNy_peaks.narrowPeak > macs2/BRD4_ChIP_IFNy_peaks.overlapped_with_IRF1_ChIP_IFNy_peaks.narrowPeak
```

　同様に、IRF1のピークのうち、BRD4のピークと重なるものを抽出する。
```
$ cd ~/chipseq
$ bedtools intersect -u -a macs2/IRF1_ChIP_IFNy_peaks.narrowPeak -b macs2/BRD4_ChIP_IFNy_peaks.narrowPeak > macs2/IRF1_ChIP_IFNy_peaks.overlapped_with_BRD4_ChIP_IFNy_peaks.narrowPeak
```

　以下のコマンドで、-a で指定したピークのうち、-bで指定したピークと重ならないものを抽出する。
```
$ cd ~/chipseq
$ bedtools intersect -v -a macs2/BRD4_ChIP_IFNy_peaks.narrowPeak -b macs2/IRF1_ChIP_IFNy_peaks.narrowPeak > macs2/BRD4_ChIP_IFNy_peaks.not_overlapped_with_IRF1_ChIP_IFNy_peaks.narrowPeak
$ bedtools intersect -v -a macs2/IRF1_ChIP_IFNy_peaks.narrowPeak -b macs2/BRD4_ChIP_IFNy_peaks.narrowPeak > macs2/IRF1_ChIP_IFNy_peaks.not_overlapped_with_BRD4_ChIP_IFNy_peaks.narrowPeak
```

　以下のコマンドで、それぞれの行数（ピーク数）を調べる。
```
$ cd ~/chipseq
$ wc -l macs2/*overlapped*.narrowPeak
    8441 macs2/BRD4_ChIP_IFNy_peaks.not_overlapped_with_IRF1_ChIP_IFNy_peaks.narrowPeak
     907 macs2/BRD4_ChIP_IFNy_peaks.overlapped_with_IRF1_ChIP_IFNy_peaks.narrowPeak
    2957 macs2/IRF1_ChIP_IFNy_peaks.not_overlapped_with_BRD4_ChIP_IFNy_peaks.narrowPeak
     909 macs2/IRF1_ChIP_IFNy_peaks.overlapped_with_BRD4_ChIP_IFNy_peaks.narrowPeak
```

### BAMファイルのBigWigファイルへの変換

```
$ cd ~/chipseq
$ mkdir deeptools    # deepToolsの出力結果を保存するディレクトリを作成する
```

```
$ cd ~/chipseq
$ bamCoverage -b bowtie2/BRD4_ChIP_IFNy.trim.uniq.bam -o deeptools/BRD4_ChIP_IFNy.trim.uniq.bw -of bigwig --normalizeUsing CPM
```

```
$ cd ~/chipseq
$ bamCoverage -b bowtie2/IRF1_ChIP_IFNy.trim.uniq.bam -o deeptools/IRF1_ChIP_IFNy.trim.uniq.bw -of bigwig --normalizeUsing CPM
$ bamCoverage -b bowtie2/Input_DNA.trim.uniq.bam -o deeptools/Input_DNA.trim.uniq.bw -of bigwig --normalizeUsing CPM
```
　以下のファイルができたことを確認する。

```
deeptools/BRD4_ChIP_IFNy.trim.uniq.bw
deeptools/IRF1_ChIP_IFNy.trim.uniq.bw
deeptools/Input_DNA.trim.uniq.bw
```

### IGV によるピークの確認（目視）

　以下のコマンドで、IGVを起動する。
```
$ igv
```

### HOMERによるモチーフ検索

　まず、HOMERの結果を保存するディレクトリを作成する。
```
$ cd ~/chipseq
$ mkdir homer
```
　以下のコマンドで、macs2/BRD4_ChIP_IFNy_summits.bedに対してHOMERが実行される。
```
$ cd ~/chipseq
$ mkdir homer/BRD4_ChIP_IFNy
$ findMotifsGenome.pl macs2/BRD4_ChIP_IFNy_summits.bed mm10 homer/BRD4_ChIP_IFNy -size 200 -mask
```

`mm10` はマウスゲノムを表す。ヒトの場合は`hg38` などにする。
`homer/` 以下に `homerResults.html` と `knownResults.html` という２つのレポートが出力される。`homerResults.html`は新規にモチーフを探索した結果を、`knownResults.html`は既知のモチーフの有無をスキャンした結果をそれぞれ記録している。

　同様に、macs2/IRF1_ChIP_IFNy_summits.bedについてもHOMERを実行する。
```
$ cd ~/chipseq
$ mkdir homer/IRF1_ChIP_IFNy
$ findMotifsGenome.pl macs2/IRF1_ChIP_IFNy_summits.bed mm10 homer/IRF1_ChIP_IFNy -size 200 -mask
```

　homer/BRD4_ChIP_IFNy および homer/IRF1_ChIP_IFNyのそれぞれに、以下のようなファイルが出力される。このうち “homerResults.html”と”knownResults.html”には、それぞれ新規モチーフの探索の結果および既知モチーフのスキャンの結果が要約されている。
```
homerMotifs.all.motifs
homerMotifs.motifs10
homerMotifs.motifs12
homerMotifs.motifs8
homerResults
homerResults.html
knownResults
knownResults.html
knownResults.txt
motifFindingParameters.txt
seq.autonorm.tsv
```


### deepToolsによるChIP-seqのリードのシグナルの分布の可視化

　以下のコマンドで、マウスの遺伝子モデルのGTFファイル(~/gencode/gencode.vM20.annotation.gtf) をregions、BRD4 ChIP-seqのbigWigファイルをscoreFile (deeptools/BRD4_ChIP_IFNy.trim.uniq.bw ) として、matrix ファイル（deepTools独自の形式のファイル）を作成する。computeMatrixコマンドの”scale-regions”というモードを使用する。

```
$ cd ~/chipseq
$ computeMatrix scale-regions \
	--regionsFileName ~/gencode/gencode.vM20.annotation.gtf \
	--scoreFileName deeptools/BRD4_ChIP_IFNy.trim.uniq.bw \
	--outFileName deeptools/BRD4_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
	--upstream 1000 --downstream 1000 \
	--skipZeros
```

　以下のコマンドで、Metagene plot を作成する。
```
$ cd ~/chipseq
$  plotProfile -m deeptools/BRD4_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
              -out deeptools/metagene_BRD4_ChIP_IFNy_gencode_vM20_gene.pdf \
              --plotTitle "GENCODE vM20 genes"
```

　次のコマンドで、遺伝子領域集合に対するChIP-seqのリードのシグナルのヒートマップを作成できる。

```
$ cd ~/chipseq
$ plotHeatmap -m deeptools/BRD4_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
              -out deeptools/heatmap_BRD4_ChIP_IFNy_gencode_vM20_gene.pdf \
              --plotTitle "GENCODE vM20 genes"
```

同様に、以下のコマンドで、IRF1のChIP-seqデータについても、metagene plotとヒートマップを作成する。

```
$ cd ~/chipseq
$ computeMatrix scale-regions \
	--regionsFileName ~/gencode/gencode.vM20.annotation.gtf \
	--scoreFileName deeptools/IRF1_ChIP_IFNy.trim.uniq.bw \
	--outFileName deeptools/IRF1_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
	--upstream 1000 --downstream 1000 \
	--skipZeros
$ plotProfile -m deeptools/IRF1_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
              -out deeptools/metagene_IRF1_ChIP_IFNy_gencode_vM20_gene.pdf \
              --plotTitle "GENCODE vM20 genes"
$ plotHeatmap -m deeptools/IRF1_ChIP_IFNy.trim.uniq.matrix_gencode_vM20_gene.txt.gz \
              -out deeptools/heatmap_IRF1_ChIP_IFNy_gencode_vM20_gene.pdf \
              --plotTitle "GENCODE vM20 genes"
```


　次に以下のコマンドで、BRD4 ChIP-seqのリードの分布 (deeptools/BRD4_ChIP_IFNy.trim.uniq.bw)がIRF1 ChIP-seqのピーク (macs2/IRF1_ChIP_IFNy_summits.bed) を中心としたときにゲノム全体としてどうなっているかをaggregation plotとして描くための準備をする。computeMatrixコマンドの”reference-point”というモードを使用する。

```
$ cd ~/
$ computeMatrix reference-point \
	--regionsFileName macs2/IRF1_ChIP_IFNy_summits.bed \
	--scoreFileName deeptools/BRD4_ChIP_IFNy.trim.uniq.bw \
	--referencePoint center \
	--upstream 1000 \
	--downstream 1000 \
	--outFileName deeptools/BRD4_ChIP_IFNy.trim.IRF1_ChIP_IFNy_summits.matrix.txt.gz \
	--skipZeros
```


　次に、以下のコマンドで、aggregation plotを作成する。

```
$ plotProfile -m deeptools/BRD4_ChIP_IFNy.trim.IRF1_ChIP_IFNy_summits.matrix.txt.gz \
              -out deeptools/aggregation_BRD4_ChIP_IFNy.trim.IRF1_ChIP_IFNy_summits.pdf \
              --regionsLabel "IRF1_ChIP_IFNy Peaks"
```

　結果は下の図のようになる。この図から、ピークの前後数百塩基の範囲にリードが集中していることがわかる。


```
$ plotHeatmap -m deeptools/BRD4_ChIP_IFNy.trim.IRF1_ChIP_IFNy_summits.matrix.txt.gz \
              -out deeptools/heatmap_BRD4_ChIP_IFNy.trim.IRF1_ChIP_IFNy_summits.pdf \
              --samplesLabel "BRD4_ChIP_IFNy" \
              --regionsLabel "IRF1_ChIP_IFNy Peaks"
```

　結果は下の図のようになる。この図から、IRF1のピークの一部について、その周辺にBRD4のリードが集中していることがわかる。


　さらに、以下のように--scoreFileName に複数のbigWigファイルを指定することで、複数のChIP-seqデータにおけるリードの分布を同時に可視化することができる。plotHeatmapで--kmeans でクラスタ数を指定することで、k-meansアルゴリズムでゲノム領域をクラスタリングして表示させることもできる。

```
$ computeMatrix scale-regions \
	--regionsFileName ../data/gencode/gencode.vM20.annotation.gtf \
	--scoreFileName deeptools/BRD4_ChIP_IFNy.trim.uniq.bw \
	deeptools/IRF1_ChIP_IFNy.trim.uniq.bw \
	--outFileName deeptools/chipseq_matrix_gencode_vM20_gene.txt.gz \
	--upstream 1000 --downstream 1000 \
	--skipZeros

$ plotHeatmap -m deeptools/chipseq_matrix_gencode_vM20_gene.txt.gz \
              -out deeptools/heatmap_BRD4_ChIP_IFNy_gencode_vM20_gene.k3.pdf \
              --kmeans 3 \
              --plotTitle "GENCODE vM20 genes"
```



### GREATによるピーク領域に対するオントロジー・パスウェイ解析

　まず 、以下のコマンドで、GREATで受け付けてもらえるように加工する。具体的には、１列目から６列目だけを抽出してBEDフォーマットにする。
```
$ cd ~/chipseq
$ cut -f 1,2,3,4,5,6 macs2/BRD4_ChIP_IFNy_peaks.narrowPeak > macs2/BRD4_ChIP_IFNy_peaks.narrowPeak.bed
$ cut -f 1,2,3,4,5,6 macs2/IRF1_ChIP_IFNy_peaks.narrowPeak > macs2/IRF1_ChIP_IFNy_peaks.narrowPeak.bed
```

　ここで、headコマンドを使いBEDファイルの先頭 10行を見て、先のコマンドの結果を確認する。
```
$ cd ~/chipseq
$ head macs2/*.narrowPeak.bed
==> macs2/BRD4_ChIP_IFNy_peaks.narrowPeak.bed <==
chr1	4807514	4808176	BRD4_ChIP_IFNy_peak_1	117	.
chr1	4857437	4857680	BRD4_ChIP_IFNy_peak_2	35	.
chr1	4857758	4858397	BRD4_ChIP_IFNy_peak_3	216	.
chr1	5018884	5019146	BRD4_ChIP_IFNy_peak_4	48	.
chr1	5019310	5019671	BRD4_ChIP_IFNy_peak_5	60	.
chr1	5022794	5023366	BRD4_ChIP_IFNy_peak_6	117	.
chr1	5082960	5083202	BRD4_ChIP_IFNy_peak_7	80	.
chr1	6214644	6215164	BRD4_ChIP_IFNy_peak_8	95	.
chr1	7088391	7088703	BRD4_ChIP_IFNy_peak_9	73	.
chr1	9747880	9748331	BRD4_ChIP_IFNy_peak_10	78	.

==> macs2/IRF1_ChIP_IFNy_peaks.narrowPeak.bed <==
chr1	3405415	3405606	IRF1_ChIP_IFNy_peak_1	95	.
chr1	3408231	3408677	IRF1_ChIP_IFNy_peak_2	586	.
chr1	6406537	6406748	IRF1_ChIP_IFNy_peak_3	143	.
chr1	6717606	6717857	IRF1_ChIP_IFNy_peak_4	207	.
chr1	7139854	7140166	IRF1_ChIP_IFNy_peak_5	25	.
chr1	7660520	7660678	IRF1_ChIP_IFNy_peak_6	107	.
chr1	9129013	9129176	IRF1_ChIP_IFNy_peak_7	95	.
chr1	9703961	9704130	IRF1_ChIP_IFNy_peak_8	90	.
chr1	9943944	9944140	IRF1_ChIP_IFNy_peak_9	35	.
chr1	10220210	10220391	IRF1_ChIP_IFNy_peak_10	88	.
```

- GREATのウェブサイト: http://great.stanford.edu

　まず、"Species Assembly"でゲノムのバージョンを選ぶ。ここでは”Mouse: NCBI build 38 (UCSC mm10, Dec/2011)”を選択する。次に、"Test regions" のBED fileには [ファイルを選択] (Choose File) をボタンをクリックして、作成したBEFファイル(IRF1_ChIP_IFNy_peaks.narrowPeak.bed) を選ぶ。


　次に、Association rule settings の [Show settings]ボタンを押すと、下図の項目が現れる。　GREATツールは、与えたBEDファイルの領域（今回はピーク領域）に対し、遺伝子を割り当てるが、その割り当て方を選択できる。今回は、[Basal plus extension]を選択し、[Submit]ボタンを押す。
　なお、エンハンサーも含まれるピーク領域を対象にする場合は、[Two nearest genes]を選択するとよい。

　[Job Description]の中のAssociated genomic regions の項目では、どのピークが何の遺伝子に割り当てられたかを示すリンク [View all genomic region-gene associations] がある。これを保存しておくと、後でどんな遺伝子が割り当てられたを確認する際に便利である。

> 2019年2月現在でヒトゲノムについてはhg19のみが使える。そのため、hg38など別のバージョンのヒトリファレンスゲノムで解析した結果でGREATを使いたい場合、
liftover (https://genome.ucsc.edu/cgi-bin/hgLiftOver)などのツールを用いてゲノムの座標をhg19に合わせて変換する必要がある。


###ChIPpeakAnnoによるピークへのアノテーション

　まず、以下のコマンドで、RStudioを起動する。

```
$ cd ~/chipseq
$ open -a RStudio
```
　なお、以下では適宜Rのパッケージをインストールする。その際には、"Update all/some/none? [a/s/n]:" と表示されたら "a" と入力してエンターを押して次に進む。

　以下では、「BRD4 ChIP-seqのピークの集合」と「IRF1 ChIP-seqのピークの集合」がどのくらい重なるかを調べる。
```
> library(ChIPpeakAnno)
> gr1 <- toGRanges("macs2/BRD4_ChIP_IFNy_peaks.narrowPeak", format="narrowPeak", header=FALSE)
> gr2 <- toGRanges("macs2/IRF1_ChIP_IFNy_peaks.narrowPeak", format="narrowPeak", header=FALSE)

> ol <- findOverlapsOfPeaks(gr1, gr2) # ピーク同士の重なりを調査する

> makeVennDiagram(ol, NameOfPeaks=c(“BRD4”, “IRF1”)) # 重なりをベン図として可視化する
```



　次に、ピークを最も近い転写開始点の遺伝子へ割り当てる。
```
> BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")    # マウスゲノムmm10の遺伝子モデルのパッケージをダウンロードする
> library(TxDb.Mmusculus.UCSC.mm10.ensGene)    # マウスゲノムmm10の遺伝子モデルのパッケージをロードする

> annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm10.ensGene)
> seqlevelsStyle(gr1) <- seqlevelsStyle(annoData)    # 染色体名のスタイルを揃える

> anno1 <- annotatePeakInBatch(gr1, AnnotationData=annoData) # ピークを最も近い転写開始点（TSS）に割り当てる
> pie1(table(anno1$insideFeature), main="BRD4")    # ピークが遺伝子からみてどの領域に置いたのかを、円グラフとして表示する
```

　上の図から、転写開始点に重なるピーク（overlapStart）が多いことがわかる。
```
> seqlevelsStyle(gr2) <- seqlevelsStyle(annoData)    # 染色体名のスタイルを揃える
> anno2 <- annotatePeakInBatch(gr2, AnnotationData=annoData) # ピークを最も近い転写開始点（TSS）に割り当てる
> pie1(table(anno2$insideFeature), main="IRF2")    # ピークが遺伝子からみてどの領域に置いたのかを、円グラフとして表示する
```

　上の図から、IRF2のピークは転写開始点の上流に多いことがわかる。

　また、以下のコマンドで、BRD4とIRF2のChIP-seqピークが重なる領域はどのような特徴を持つかを調べることができる。
```
> overlaps <- ol$peaklist[["gr1///gr2"]]
> aCR <- assignChromosomeRegion(overlaps, nucleotideLevel=FALSE,
                           precedence=c("Promoters", "immediateDownstream",
                                         "fiveUTRs", "threeUTRs",
                                         "Exons", "Introns"),
                           TxDb=TxDb.Mmusculus.UCSC.mm10.ensGene)
> pie1(aCR$percentage, main="BRD4 & IRF1")
```

　上のアノテーションには遺伝子IDは含まれているが、遺伝子名が含まれていない。遺伝子名も対応づいていたほうが便利なことが多いため、遺伝子IDに対応する遺伝子名を追加する。
```
> BiocManager::install("EnsDb.Mmusculus.v79")
> library(EnsDb.Mmusculus.v79)

> anno1$feature[is.na(anno1$feature)] <- "."  # エラーを避けるために NA をピリオドに変える
> anno1$geneName <- mapIds(EnsDb.Mmusculus.v79, keys=anno1$feature, column = "GENENAME", keytype="GENEID")

> anno1[1:2]
GRanges object with 2 ranges and 14 metadata columns:
                                           seqnames          ranges strand |
                                              <Rle>       <IRanges>  <Rle> |
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903     chr1 4807514-4808176      * |
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813     chr1 4857437-4857680      * |
                                               score signalValue    pValue
                                           <integer>   <numeric> <numeric>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903       117     8.68079  15.22217
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813        35     4.76915   6.42381
                                              qValue                  peak
                                           <numeric>           <character>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903   11.7561 BRD4_ChIP_IFNy_peak_1
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813   3.56659 BRD4_ChIP_IFNy_peak_2
                                                      feature start_position
                                                  <character>      <integer>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903 ENSMUSG00000025903        4807788
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813 ENSMUSG00000033813        4857814
                                           end_position feature_strand
                                              <integer>    <character>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903      4886770              +
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813      4897909              +
                                           insideFeature distancetoFeature
                                                <factor>         <numeric>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903  overlapStart              -274
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813      upstream              -377
                                           shortestDistance
                                                  <integer>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903              274
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813              134
                                           fromOverlappingOrNearest    geneName
                                                        <character> <character>
  BRD4_ChIP_IFNy_peak_1.ENSMUSG00000025903          NearestLocation      Lypla1
  BRD4_ChIP_IFNy_peak_2.ENSMUSG00000033813          NearestLocation       Tcea1
  -------
  seqinfo: 23 sequences from an unspecified genome; no seqlengths
```

　以下のコマンドで、結果をタブ区切りファイルとして保存する。
```
if(!dir.exists("ChIPpeakAnno")) dir.create("ChIPpeakAnno")
df_anno1 <- as.data.frame(anno1)
write.table(df_anno1, "ChIPpeakAnno/BRD4_ChIP_IFNy_peaks.annot.txt", sep="\t", quote=F)
```
