SELECT CONCAT('aws s3 cp --no-sign-request s3://1000genomes-dragen-3.7.6/data/individuals/hg38-graph-based', relative_path) FROM DIRECTORY( @dragen_all )
where --relative_path rlike '/HG0011[0-7].*.hard-filtered.vcf.gz'
relative_path rlike '/*.*.hard-filtered.vcf.gz'
order by SIZE desc
limit 100 ;