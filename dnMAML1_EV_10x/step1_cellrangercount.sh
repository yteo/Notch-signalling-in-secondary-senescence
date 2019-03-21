# used cellranger-2.0.1
cellranger count --id=dnMAML1_mVenus \
--transcriptome=/data/hg19_mVenus_puro_neo \
--fastqs=/data/dMAML1 \
--sample=1_11376PLpool10__dMAM,2_11376PLpool10__dMAM \
--localmem=80 \
--localcores=20


cellranger count --id=EV_mVenus \
--transcriptome=/data/hg19_mVenus_puro_neo \
--fastqs=/data/EV \
--sample=1_11376PLpool10__EV,2_11376PLpool10__EV \
--localmem=80 \
--localcores=20