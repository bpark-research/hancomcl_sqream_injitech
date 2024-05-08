# hancomcl_sqream_injitech

(1) hail로 1KG 100개 vcf 합치는 작업: [hail/bpark]</br>
(2) sqream python 연결: [sqream/bpark]</br>

문서작성: https://docs.google.com/document/d/1dxSrIM52xJ5PlqMigxYxRSNkl1ufW0I1-BFfjq4M9nE/edit?usp=sharing

### Reference
vcf file에 대한 문서 [https://samtools.github.io/hts-specs/VCFv4.2.pdf]

python vcf 라이브러리</br>
[https://vcfpy.readthedocs.io/en/stable/]</br>
[https://scikit-allel.readthedocs.io/en/stable/index.html] 개발 중단 -> sgkit</br>
[https://github.com/sgkit-dev/sgkit]</br>

scalable gVCF merging and joint variant calling for population sequencing projects [https://github.com/dnanexus-rnd/GLnexus]

레드헷에 bpark 계정</br>
mkdir -p ~/miniconda3</br>
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh</br>
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3</br>
rm -rf ~/miniconda3/miniconda.sh</br>
conda -n hail python=3.9</br>
conda activate hail</br>

[hail/bpark]: hail/bpark.ipynb
[sqream/bpark]: sqream/bpark.ipynb
