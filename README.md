# hancomcl_sqream_injitech

(1) hail로 3개 vcf 합치는 작업: [hail/bpark]
(2) sqream python 연결: [sqream/bpark]

vcf file에 대한 문서 [https://samtools.github.io/hts-specs/VCFv4.2.pdf]

python vcf 라이브러리
[https://vcfpy.readthedocs.io/en/stable/]
[https://scikit-allel.readthedocs.io/en/stable/index.html] 개발 중단 -> sgkit
[https://github.com/sgkit-dev/sgkit]

scalable gVCF merging and joint variant calling for population sequencing projects [https://github.com/dnanexus-rnd/GLnexus]

레드헷에 bpark 계정
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

[hail/bpark]: hail/bpark.ipynb
[sqream/bpark]: sqream/bpark.ipynb